/*!
* \file main.cc
* \brief Main file of the TEST_XTRX_OCTO program.
* \author Javier Arribas, 2021. jarribas(at)cttc.es
*
* It sets up the logging system, creates a the main thread,
* makes it run, and releases memory back when the main thread has ended.
*
* -------------------------------------------------------------------------
*
* Copyright (C) 2021 (see AUTHORS file for a list of contributors)
*
* This file is part of TEST_XTRX_OCTO.
*
*/

#ifndef TEST_XTRX_OCTO_VERSION
#define TEST_XTRX_OCTO_VERSION "0.0.1"
#endif

//#include "concurrent_queue.h"
#include "test_xtrx_octo_flags.h"
#include <xtrx_api.h>
#include "xtrx_octo_api.h"
#include <signal.h>
#include <semaphore.h>
#include <boost/exception/diagnostic_information.hpp>  // for diagnostic_information
#include <boost/exception/exception.hpp>               // for exception
#include <boost/thread/exceptions.hpp>                 // for thread_resource_error
#include <gflags/gflags.h>
#include <gflags/gflags_declare.h>
#include <glog/logging.h>  // for FLAGS_log_dir
#include <chrono>
#include <exception>  // for exception
#include <iostream>
#include <sys/stat.h>
#include <thread>



enum {
	BUF_TO_FLUSH = 32,
	MAX_DEVS     = 8,
};

struct stream_data {
	struct xtrx_dev *dev;

	uint64_t cycles;
	uint64_t samples_per_cyc;
	double ramplerate;
	unsigned slice_sz;

	bool mux_demux; // Data multiplexing/demultiplexing
	bool siso;      // Device operatig mode
	bool flush_data;
	unsigned sample_host_size;

	unsigned out_overruns;

	size_t out_tm_diff;

	uint64_t out_samples_per_dev;
	uint64_t out_time;
};

struct rx_flush_data {
	void *buffer_ptr[MAX_DEVS*2][BUF_TO_FLUSH];
	FILE *out_files[MAX_DEVS*2];

	unsigned num_chans;
	unsigned buf_ptr;
	unsigned buf_max;

	unsigned wr_sz;

	unsigned put_cnt;
	unsigned get_cnt;

	bool stop;
	bool fiop;

	sem_t sem_read_rx;
	sem_t sem_read_wr;
};

static uint64_t grtime(void)
{
	struct timespec ts;
	clock_gettime(CLOCK_MONOTONIC_RAW, &ts);

	return (uint64_t)ts.tv_sec * (uint64_t)1000000000 + (uint64_t)ts.tv_nsec;
}


char* alloc_flush_data_bufs(struct rx_flush_data *pd,
							unsigned buffer_count,
							unsigned num_chans,
							unsigned slice_sz,
							const char* tag,
							const char* basename,
							const char* fattr,
							bool rd)
{
	pd->buf_max = buffer_count;
	pd->buf_ptr = 0;
	pd->num_chans = num_chans;
	pd->stop = false;
	pd->fiop = basename != NULL;
	pd->put_cnt = 0;
	pd->get_cnt = 0;
	pd->wr_sz = slice_sz;
	sem_init(&pd->sem_read_rx, 0, (rd) ? (buffer_count) : 0);
	sem_init(&pd->sem_read_wr, 0, (!rd) ? (buffer_count) : 0);

	size_t rx_bufsize  = num_chans * pd->buf_max * pd->wr_sz;
	size_t rx_bufslice = pd->wr_sz;

	char* out_buffs = (char*)malloc(rx_bufsize);
	if (out_buffs == NULL) {
		fprintf(stderr, "Unable to create %s buffers (size=%.3fMB)\n", tag,
				(float)rx_bufsize / 1024 / 1024);
		return NULL;
	}

	char* bptr = out_buffs;
	for (unsigned p = 0; p < num_chans; p++) {
		char filename_buf[256];
		const char* filename = (num_chans == 1) ? basename : filename_buf;
		if (num_chans > 1) {
			snprintf(filename_buf, sizeof(filename_buf), "%s_%d",
					 basename, p);
		}
		if (basename) {
			pd->out_files[p] = fopen(filename, fattr);
			if (pd->out_files[p] == NULL) {
				fprintf(stderr, "Unable to open file: %s! error: %d\n",
						filename, errno);
				exit(EXIT_FAILURE);
			}
		}

		for (unsigned b = 0; b < pd->buf_max; b++) {
			pd->buffer_ptr[p][b] = bptr;
			bptr += rx_bufslice;
		}
	}

	return out_buffs;
}

int stream_rx(struct stream_data* sdata,
			  struct rx_flush_data* rxd)
{
	int res = 0;
	void* stream_buffers[2 * MAX_DEVS];
	unsigned buf_cnt = rxd->num_chans;
	uint64_t zero_inserted = 0;
	int overruns = 0;
	uint64_t rx_processed = 0;
	uint64_t sp = grtime();

	uint64_t st = sp;
	uint64_t abpkt = sp;
	unsigned binx = 0;
	unsigned dev_count = rxd->num_chans / (sdata->mux_demux ? 2 : 1);

	fprintf(stderr, "RX CYCLES=%" PRIu64 " SAMPLES=%" PRIu64 " SLICE=%u (PARTS=%" PRIu64 ")\n",
			sdata->cycles,
			sdata->samples_per_cyc, sdata->slice_sz,
			sdata->samples_per_cyc / sdata->slice_sz);


	for (uint64_t p = 0; p < sdata->cycles; p++) {
		// Get new buffer to writing to
		if (sdata->flush_data) {
			res = sem_trywait(&rxd->sem_read_rx);
			if (res) {
				fprintf(stderr, "RX rate is too much; RX_TO_FILE thread is lagging behind\n");
				goto falied_stop_rx;
			}
		}

		for (unsigned bc = 0; bc < buf_cnt; bc++) {
			stream_buffers[bc] = rxd->buffer_ptr[bc][binx];
		}

		for (uint64_t h = 0; h < sdata->samples_per_cyc / sdata->slice_sz; h++) {
			xtrx_recv_ex_info_t ri;
			size_t rem = sdata->samples_per_cyc - h * sdata->slice_sz;
			if (rem > sdata->slice_sz)
				rem = sdata->slice_sz;

			ri.samples = rem / ((sdata->mux_demux) ? 2 : 1);
			ri.buffer_count = buf_cnt;
			ri.buffers = stream_buffers;
			ri.flags = 0;
			int rx_rep_gtime=0;
			if (rx_rep_gtime) {
				ri.flags |= RCVEX_REPORT_GTIME;
			}

			uint64_t sa = grtime();
			uint64_t da = sa - sp;
			res = xtrx_recv_sync_ex(sdata->dev, &ri);
			sp = grtime();
			uint64_t sb = sp - sa;

			rx_processed += ri.out_samples * ((sdata->mux_demux) ? 2 : 1);

			abpkt += 1e9 * ri.samples * ri.buffer_count / dev_count / sdata->ramplerate / (sdata->siso ? 1 : 2);
			unsigned s_logp = 1;
			if (s_logp == 0 || p % s_logp == 0)
				fprintf(stderr, "PROCESSED RX SLICE %" PRIu64 " /%" PRIu64 ":"
								" res %d TS:%8" PRIu64 " %c%c  %6" PRId64 " us"
								" DELTA %6" PRId64 " us LATE %6" PRId64 " us"
								" %d samples\n",
						p, h, res, ri.out_first_sample,
						(ri.out_events & RCVEX_EVENT_OVERFLOW)    ? 'O' : ' ',
						(ri.out_events & RCVEX_EVENT_FILLED_ZERO) ? 'Z' : ' ',
						sb / 1000, da / 1000, (int64_t)(sp - abpkt) / 1000, ri.out_samples);
			if (res) {
				fprintf(stderr, "Failed xtrx_recv_sync: %d\n", res);
				goto falied_stop_rx;
			}

			if (ri.out_events & RCVEX_EVENT_OVERFLOW) {
				overruns++;
				zero_inserted += (ri.out_resumed_at - ri.out_overrun_at);
			}
			for (unsigned bc = 0; bc < buf_cnt; bc++) {
				stream_buffers[bc] += ri.samples * sdata->sample_host_size;
			}
		}

		if (sdata->flush_data) {
			binx = (binx + 1) % rxd->buf_max;
			rxd->put_cnt++;

			res = sem_post(&rxd->sem_read_wr);
			if (res) {
				fprintf(stderr, "RX unable to post buffers!\n");
				goto falied_stop_rx;
			}
		}
	}
falied_stop_rx:
	sdata->out_overruns = overruns;
	sdata->out_samples_per_dev = rx_processed;
	sdata->out_tm_diff = grtime() - st;
	return res;
}

void* thread_rx_to_file(void* obj)
{
	struct rx_flush_data *rxd = (struct rx_flush_data *)obj;
	int res;

	fprintf(stderr, "RX_TO_FILE: thread start\n");

	for (;;) {
		res = sem_wait(&rxd->sem_read_wr);
		if (res) {
			fprintf(stderr, "RX_TO_FILE: sem wait error!\n");
			return 0;
		}
		if ((rxd->put_cnt == rxd->get_cnt) && rxd->stop) {
			fprintf(stderr, "RX_TO_FILE: thread exit: %d:%d\n", rxd->put_cnt, rxd->get_cnt);
			return 0;
		}

		for (unsigned i = 0; i <rxd->num_chans; i++) {
			size_t ret = fwrite(rxd->buffer_ptr[i][rxd->buf_ptr], rxd->wr_sz, 1, rxd->out_files[i]);
			if (ret != 1) {
				fprintf(stderr, "RX_TO_FILE: write error %u != expected %u\n", (unsigned)ret, rxd->wr_sz);
			}
		}
		sem_post(&rxd->sem_read_rx);
		rxd->buf_ptr = (rxd->buf_ptr + 1) % rxd->buf_max;

		//fprintf(stderr, "RX_TO_FILE: buffer %u written\n", rxd->get_cnt);
		rxd->get_cnt++;
	}
	return 0;
}

bool file_exist(const std::string& name)
{
    return (access(name.c_str(), F_OK) != -1);
}

int main(int argc, char** argv)
{
    const std::string intro_help(
        std::string("\nTEST_XTRX_OCTO_VERSION is a command line tool to manage the ESA MultiRF GNSS Dongle and capture baseband samples with uBlox timing information\n") +
        "Copyright (C) 2021 (see AUTHORS file for a list of contributors)\n" +
        "This program comes with ABSOLUTELY NO WARRANTY;\n" +
        "See COPYING file to see a copy of the General Public License\n \n");

    const std::string version(TEST_XTRX_OCTO_VERSION);
#ifdef GFLAGS_NAMESPACE
    GFLAGS_NAMESPACE::SetUsageMessage(intro_help);
    GFLAGS_NAMESPACE::SetVersionString(version);
    GFLAGS_NAMESPACE::ParseCommandLineFlags(&argc, &argv, true);
#else
    google::SetUsageMessage(intro_help);
    google::SetVersionString(gnss_sdr_version);
    google::ParseCommandLineFlags(&argc, &argv, true);
#endif

    std::cout << "Initializing DONGLE-CAPTURE v" << version << " ... Please wait." << std::endl;
    FLAGS_log_dir = "./";
    google::InitGoogleLogging(argv[0]);
    std::cout << "Logging will be written at " << FLAGS_log_dir << '\n';

    //find XTRX devices
	const unsigned MAX_DEVS = 32;
	xtrx_device_info_t devs[MAX_DEVS];
	int res = xtrx_discovery(devs, MAX_DEVS);
	xtrx_dev* dev=NULL;
	if (res > 0)
	{
		std::vector<std::string> xtrx_devs_str;
		for (int i = 0; i < res; i++) {
			xtrx_devs_str.push_back(devs[i].uniqname);
			std::cout<<"["<<i<<"] found device: "<<devs[i].uniqname<<"\n";

			}

		//int res_open = xtrx_open_string(xtrx_devs_str[0].c_str(), &dev);
		//note: in XYNC octo card, there is one XTRX that have access to the common front-end logic.
		//in my setup, this is the xtrx4 the one that have access...
		int res_open = xtrx_open_string("/dev/xtrx4;"
				"/dev/xtrx0;/dev/xtrx1;/dev/xtrx2;/dev/xtrx3;/dev/xtrx5;/dev/xtrx6;/dev/xtrx7"
				";;fe=octoCAL;loglevel=4", &dev);
		//int res_open = xtrx_open_string("/dev/xtrx4;;fe=octoCAL;loglevel=2", &dev);
		//int res_open = xtrx_open_string(xtrx_devs_str[0].c_str(), &dev);
		//int res_open = xtrx_open_string(";;loglevel=4", &dev);
		int dev_count=8;
		//int res_open = xtrx_open_string("/dev/xtrx4;;fe=octoRFX6:clk;loglevel=4", &dev);
		//int xtrx_open_multi(const xtrx_open_multi_info_t *dinfo, struct xtrx_dev** outdev)

		if (res_open >= 0)
		{
			int n_subdev = res;
			std::cout<<"Device "<<xtrx_devs_str[0]<<" openned with "<<n_subdev<<" subdevices\n";

			//OCTO PCI CARD Calibration
			//int tmpres  = xtrx_set_ref_clk(dev, 0, XTRX_CLKSRC_EXT_W1PPS_SYNC);
			int tmpres  = xtrx_set_ref_clk(dev, 0, XTRX_CLKSRC_INT);
			std::cout <<" xtrx_set_ref_clk returned "<<tmpres<<"\n";
			if (res_open >= 0)
			{
				std::cout<<"External clock source and synchronization signal set ok\n";

				//configure stream
				xtrx_channel_t ch = XTRX_CH_ALL;
				xtrx_wire_format_t rx_wire_fmt = XTRX_WF_16;
				xtrx_host_format_t rx_host_fmt = XTRX_IQ_INT16;
				xtrx_wire_format_t tx_wire_fmt = XTRX_WF_16;
				xtrx_host_format_t tx_host_fmt = XTRX_IQ_INT16;

				//
				// Set XTRX parameters
				//
				double rxsamplerate = 4.0e6, actual_rxsample_rate = 0;
				double txsamplerate = 0, actual_txsample_rate = 0;
				double rxfreq = 900e6, rxactualfreq = 0;
				double txfreq = 450e6, txactualfreq = 0;
				double rxbandwidth = 2e6, actual_rxbandwidth = 0;
				int samples_flag=0;
				double master;
				double master_in=0;

				//fprintf(stderr,"%f,%f,%f,%i ",master_in, rxsamplerate, txsamplerate, samples_flag);
				res = xtrx_set_samplerate(dev, master_in, rxsamplerate, txsamplerate, samples_flag,
										  &master, &actual_rxsample_rate, &actual_txsample_rate);
				if (res) {
					fprintf(stderr, "Failed xtrx_set_samplerate: %d\n", res);
				}

				fprintf(stderr, "Master: %.3f MHz; RX rate: %.3f MHz; TX rate: %.3f MHz\n",
						master / 1e6,
						actual_rxsample_rate / 1e6,
						actual_txsample_rate / 1e6);

//				if (vio) {
//					xtrx_val_set(dev, XTRX_TRX, XTRX_CH_ALL, XTRX_LMS7_VIO, vio);
//				}

					xtrx_set_antenna(dev, XTRX_RX_AUTO);

					res = xtrx_tune(dev, XTRX_TUNE_RX_FDD, rxfreq, &rxactualfreq);
					if (res) {
						fprintf(stderr, "Failed xtrx_tune: %d\n", res);
					}
					fprintf(stderr, "RX tunned: %f\n", rxactualfreq);

			#if 0
					if (rxfreq < 900e6) {
						xtrx_set_antenna(dev, XTRX_RX_L);
					} else if (rxfreq > 2300e6) {
						xtrx_set_antenna(dev, XTRX_RX_H);
					} else {
						xtrx_set_antenna(dev, XTRX_RX_W);
					}
			#endif

					res = xtrx_tune_rx_bandwidth(dev, ch, rxbandwidth, &actual_rxbandwidth);
					if (res) {
						fprintf(stderr, "Failed xtrx_tune_rx_bandwidth: %d\n", res);
						//goto falied_tune;
					}
					fprintf(stderr, "RX bandwidth: %f\n", actual_rxbandwidth);

					double rxgain_lna = 15;
					double actuallnagain;
					res = xtrx_set_gain(dev, ch, XTRX_RX_LNA_GAIN, rxgain_lna, &actuallnagain);
					if (res) {
						fprintf(stderr, "Failed xtrx_set_gain(LNA): %d\n", res);
					}
					fprintf(stderr, "RX LNA gain: %f\n", actuallnagain);
					int rxgain_pga = 0;
					int rxgain_tia = 9;

					res = xtrx_set_gain(dev, ch, XTRX_RX_PGA_GAIN, rxgain_pga, &actuallnagain);
					if (res) {
						fprintf(stderr, "Failed xtrx_set_gain(PGA): %d\n", res);
					}
					fprintf(stderr, "RX PGA gain: %f\n", actuallnagain);

					res = xtrx_set_gain(dev, ch, XTRX_RX_TIA_GAIN, rxgain_tia, &actuallnagain);
					if (res) {
						fprintf(stderr, "Failed xtrx_set_gain(TIA): %d\n", res);
					}
					fprintf(stderr, "RX TIA gain: %f\n", actuallnagain);

					xtrx_direction_t dir = XTRX_RX;
					//dir |= XTRX_RX;
					int loopback=0;
					int rxlfsr=0;
					int rx_tst_a = 0;
					int rx_tst_b = 0;
					int tx_tst_a = 0;
					int tx_tst_b = 0;
					int tx_siso = 0;
					int rx_siso = 0;
					int tx_swap_ab = 0;
					int rx_swap_ab = 0;
					int tx_swap_iq = 0;
					int rx_swap_iq = 0;
					int rx_packet_size = 0;
					int tx_packet_size = 0;
					unsigned rx_skip = 8192;
					unsigned tx_skip = 8192;
					int gmode = 0;

					xtrx_run_params_t params;
						params.dir = dir;
						params.nflags = (loopback) ? XTRX_RUN_DIGLOOPBACK :
													 (rxlfsr)   ? XTRX_RUN_RXLFSR : 0;
						params.rx.wfmt = rx_wire_fmt;
						params.rx.hfmt = rx_host_fmt;
						params.rx.chs = ch;
						params.rx.flags = (rx_tst_a    ? XTRX_RSP_TEST_SIGNAL_A : 0) |
								(rx_tst_b    ? XTRX_RSP_TEST_SIGNAL_B : 0) |
								(rx_siso     ? XTRX_RSP_SISO_MODE     : 0) |
								(rx_swap_ab  ? XTRX_RSP_SWAP_AB       : 0) |
								(rx_swap_iq  ? XTRX_RSP_SWAP_IQ       : 0);
						params.rx.paketsize = rx_packet_size / ((rx_siso) ? 1 : 2);
						params.tx.wfmt = tx_wire_fmt;
						params.tx.hfmt = tx_host_fmt;
						params.tx.chs = ch;
						params.tx.flags = (tx_tst_a    ? XTRX_RSP_TEST_SIGNAL_A : 0) |
								(tx_tst_b    ? XTRX_RSP_TEST_SIGNAL_B : 0) |
								(tx_siso     ? XTRX_RSP_SISO_MODE     : 0) |
								(tx_swap_ab  ? XTRX_RSP_SWAP_AB       : 0) |
								(tx_swap_iq  ? XTRX_RSP_SWAP_IQ       : 0);
						params.tx.paketsize = tx_packet_size / ((tx_siso) ? 1 : 2);
						params.rx_stream_start = rx_skip;
						params.tx_repeat_buf = NULL; //(tx_repeat_mode) ? BUF : NULL;
						params.gtime.sec = 11;
						params.gtime.nsec = 750000000; //250 ms delay
						if (gmode > 0) {
							params.nflags |= XTRX_RUN_GTIME;
						}


						xtrx_stop(dev, XTRX_TRX);
						xtrx_stop(dev, XTRX_TRX);

						res = xtrx_run_ex(dev, &params);
						if (res) {
							fprintf(stderr, "Failed xtrx_run: %d\n", res);
						}


				//create RX thread
				pthread_t rx_write_thread;
				struct rx_flush_data rxd;
				char* out_buffs = NULL;
				unsigned rx_sample_host_size = sizeof(float) * 2;
				uint64_t samples = 16384;

				int mimomode=1;

				struct stream_data sd_rx;
				sd_rx.dev = dev;
				sd_rx.cycles = 1;
				sd_rx.samples_per_cyc = samples;
				sd_rx.slice_sz = samples;
				sd_rx.ramplerate = actual_rxsample_rate;
				sd_rx.mux_demux = mimomode;
				sd_rx.siso = rx_siso;
				sd_rx.flush_data = 1;
				sd_rx.sample_host_size = rx_sample_host_size;


				int rx_buffer_count=16;
					out_buffs = alloc_flush_data_bufs(&rxd,
													  rx_buffer_count,
													  dev_count * (mimomode ? 2 : 1),
													  rx_sample_host_size * samples / (mimomode ? 2 : 1),
													  "RX",
													  "rx_signal.dat",
													  "wb",
													  true);
					if (out_buffs == NULL) {
						std::cout<<"alloc_flush_data_bufs create fail\n";
					}
					res = pthread_create(&rx_write_thread, NULL, thread_rx_to_file, &rxd);
					if (res) {
						std::cout<<"pthread create fail\n";
					}

					stream_rx(&sd_rx, &rxd);
					uint64_t rx_tm = 0;
					uint64_t rx_processed = 0;
					rx_tm = sd_rx.out_tm_diff;
					rx_processed = sd_rx.out_samples_per_dev;

					fprintf(stderr, "RX STAT Overruns:%d\n", sd_rx.out_overruns);

						rxd.stop = true;
						sem_post(&rxd.sem_read_wr);
						pthread_join(rx_write_thread, NULL);

						fprintf(stderr, "Success!\n");


						unsigned rxchs = (rx_siso ? 1 : 2);

						double rx_t = (double)rx_processed*1000/rx_tm;
						double rx_w = (double)rx_processed*(XTRX_WF_16+1)*1000/rx_tm;

						fprintf(stderr, "Processed %d devs, each: RX %d x %.3f = %.3f MSPS (WIRE: %f) \n",
								dev_count,
								rxchs, rx_t / rxchs, rx_t, rx_w);


			}
			//int tmpres = xtrx_octo_set_cal_path(dev, true);

		}else{
			std::cout<<"Unable to open XTRX device!\n";
		}

		if (dev) {
			xtrx_stop(dev, XTRX_TRX);
			xtrx_close(dev);
			dev = NULL;
		}
	}else{
		std::cout<<"No XTRX devices found on this system...\n";
	}

    google::ShutDownCommandLineFlags();
    std::cout << "DONGLE-CAPTURE program ended." << std::endl;
    int return_code = 0;

    return return_code;
}
