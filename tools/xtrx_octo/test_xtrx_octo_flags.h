#ifndef TEST_XTRX_OCTO_FLAGS_H
#define TEST_XTRX_OCTO_FLAGS_H


#include <gflags/gflags.h>
#include <cstdint>

DEFINE_uint64(freq, 1575420000, "SDR RF frequency [Hz]. Range: 100000-3800000000");

DEFINE_int32(gain, 60, "Controls combined RX gain settings. Gain range must be [0:73] [dB].");

DEFINE_int32(fs, 30720000, "LimeSDR-USB sample rate must be no more than 61.44e6 S/s [sps].");

DEFINE_int32(analogbw, 0, "Enter analog filter bandwidth for each channel. Analog filter is off if bandwidth is set to 0. Analog filter bandwidth range must be [1.5e6,130e6] Hz.");

DEFINE_int32(digitbw, 0, "Enter digital filter bandwidth for each channel. Digital filter if off if bandwidth is set to 0. Bandwidth should not be higher than sample rate.");

DEFINE_int32(ch, 1, "Number of channels to capture. [1,2]");

DEFINE_int32(ch0_antenna_port, 2,
    "Select CH0 antenna port:\n"
    "LMS_PATH_NONE = 0, <No active path (RX or TX)\n"
    "LMS_PATH_LNAH = 1, <RX LNA_H port\n"
    "LMS_PATH_LNAL = 2, <RX LNA_L port\n"
    "LMS_PATH_LNAW = 3, <RX LNA_W port\n"
    "LMS_PATH_TX1 = 1,  <TX port 1\n"
    "LMS_PATH_TX2 = 2,   <TX port 2\n"
    "LMS_PATH_AUTO = 255, <Automatically select port (if supported)");

DEFINE_int32(ch1_antenna_port, 2,
    "Select CH0 antenna port:\n"
    "LMS_PATH_NONE = 0, <No active path (RX or TX)\n"
    "LMS_PATH_LNAH = 1, <RX LNA_H port\n"
    "LMS_PATH_LNAL = 2, <RX LNA_L port\n"
    "LMS_PATH_LNAW = 3, <RX LNA_W port\n"
    "LMS_PATH_TX1 = 1,  <TX port 1\n"
    "LMS_PATH_TX2 = 2,   <TX port 2\n"
    "LMS_PATH_AUTO = 255, <Automatically select port (if supported)");

DEFINE_int32(duration, 0, "Number of seconds to capture [s]. Set 0 to capture until q key + enter is pressed");

DEFINE_string(format, "int16", "Capture format: interleaved int8H (MSByte) or int8L (LSByte) [I,Q,I,Q..] or interleaved int16 [I,Q,I,Q...]. Notice that the LimeSDR ADC resolution is 12 bits and capturing in int8 will truncate the LSB or the MSB bits, depending on the chosen option.");

DEFINE_string(capfile, "",
    "Capture filename base string. The utility will produce the following files: \n \"capfilename\"_ch0.bin \n \"capfilename\"_ch0_time.bin\n \"capfilename\"_ch1.bin\n \"capfilename\"_ch1_time.bin\n \"capfilename\"_gnss_rx1.ubx\n \"capfilename\"_gnss_rx2.ubx.\n If the file is notspecified,"
    " the utility will usethe system dateand time togenerate unique filename.");

DEFINE_bool(ppsmode, true, "Enable custom FPGA PPS mode (required for timestamping and enabled by default)");
DEFINE_bool(capture_samples, true, "Enable samples capture record to file (enabled by default)");

DEFINE_int32(filelog_level, 0, "Set log file debug verbosity level: [0,1] ");
DEFINE_int32(displaylog_level, 0, "Set display debug verbosity level: [0,1] ");

DEFINE_double(extclk_mhz, 0.0, "Enable external reference clock input and set its frequency [MHz]");

#endif  // TEST_XTRX_OCTO_FLAGS_H
