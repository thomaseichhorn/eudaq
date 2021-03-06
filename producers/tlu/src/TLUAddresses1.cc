#include "tlu/TLUAddresses.hh"
#include "tlu/TLU_address_map.h"

#define REGISTERED_BUFFER_POINTER_ADDRESS_BASE                                 \
  REGISTERED_BUFFER_POINTER_ADDRESS_0
#define REGISTERED_BUFFER_POINTER_ADDRESS_BYTES 2
#define REGISTERED_TIMESTAMP_ADDRESS_BASE REGISTERED_TIMESTAMP_ADDRESS_0
#define REGISTERED_TIMESTAMP_ADDRESS_BYTES 8
#define REGISTERED_TRIGGER_COUNTER_ADDRESS_BASE                                \
  REGISTERED_TRIGGER_COUNTER_ADDRESS_0
#define REGISTERED_TRIGGER_COUNTER_ADDRESS_BYTES 4
#define BUFFER_POINTER_ADDRESS_BASE BUFFER_POINTER_ADDRESS_0
#define BUFFER_POINTER_ADDRESS_BYTES 2
#define TIMESTAMP_ADDRESS_BASE TIMESTAMP_ADDRESS_0
#define TIMESTAMP_ADDRESS_BYTES 8
#define TRIGGER_COUNTER_ADDRESS_BASE TRIGGER_COUNTER_ADDRESS_0
#define TRIGGER_COUNTER_ADDRESS_BYTES 4
#define TRIGGER_IN0_COUNTER_BASE TRIGGER_IN0_COUNTER_0
#define TRIGGER_IN0_COUNTER_BYTES 2
#define TRIGGER_IN1_COUNTER_BASE TRIGGER_IN1_COUNTER_0
#define TRIGGER_IN1_COUNTER_BYTES 2
#define TRIGGER_IN2_COUNTER_BASE TRIGGER_IN2_COUNTER_0
#define TRIGGER_IN2_COUNTER_BYTES 2
#define TRIGGER_IN3_COUNTER_BASE TRIGGER_IN3_COUNTER_0
#define TRIGGER_IN3_COUNTER_BYTES 2
#define REGISTERED_PARTICLE_COUNTER_ADDRESS_BASE                               \
  REGISTERED_PARTICLE_COUNTER_ADDRESS_0
#define REGISTERED_PARTICLE_COUNTER_ADDRESS_BYTES 4

#define I2C_SDA_OUT_BIT EUDAQ_TLU_MISSING
#define I2C_SDA_IN_BIT EUDAQ_TLU_MISSING
#define I2C_SCL_OUT_BIT EUDAQ_TLU_MISSING
#define I2C_SCL_IN_BIT EUDAQ_TLU_MISSING
#define DUT_I2C_BUS_SELECT_ADDRESS EUDAQ_TLU_MISSING
#define DUT_I2C_BUS_DATA_ADDRESS EUDAQ_TLU_MISSING
#define I2C_BUS_MOTHERBOARD EUDAQ_TLU_MISSING
#define I2C_BUS_HDMI EUDAQ_TLU_MISSING
#define I2C_BUS_LEMO EUDAQ_TLU_MISSING
#define I2C_BUS_DISPLAY EUDAQ_TLU_MISSING
#define I2C_BUS_PMT_DAC EUDAQ_TLU_MISSING
#define I2C_BUS_MOTHERBOARD_LED_IO EUDAQ_TLU_MISSING
#define I2C_BUS_MOTHERBOARD_TRIGGER__ENABLE_IPSEL_IO EUDAQ_TLU_MISSING
#define I2C_BUS_MOTHERBOARD_RESET_ENABLE_IO EUDAQ_TLU_MISSING
#define I2C_BUS_MOTHERBOARD_FRONT_PANEL_IO EUDAQ_TLU_MISSING
#define I2C_BUS_MOTHERBOARD_LCD_IO EUDAQ_TLU_MISSING
#define I2C_BUS_LEMO_RELAY_IO EUDAQ_TLU_MISSING
#define I2C_BUS_LEMO_LED_IO EUDAQ_TLU_MISSING
#define I2C_BUS_LEMO_LED_IO_VB EUDAQ_TLU_MISSING
#define I2C_BUS_LEMO_DAC EUDAQ_TLU_MISSING
#define STROBE_PERIOD_ADDRESS_0 EUDAQ_TLU_MISSING
#define STROBE_PERIOD_ADDRESS_1 EUDAQ_TLU_MISSING
#define STROBE_PERIOD_ADDRESS_2 EUDAQ_TLU_MISSING
#define STROBE_WIDTH_ADDRESS_0 EUDAQ_TLU_MISSING
#define STROBE_WIDTH_ADDRESS_1 EUDAQ_TLU_MISSING
#define STROBE_WIDTH_ADDRESS_2 EUDAQ_TLU_MISSING
#define STROBE_ENABLE_ADDRESS EUDAQ_TLU_MISSING
#define HANDSHAKE_MODE_ADDRESS EUDAQ_TLU_MISSING
#define BUFFER_STOP_MODE_ADDRESS EUDAQ_TLU_MISSING
#define WRITE_TRIGGER_BITS_MODE_ADDRESS EUDAQ_TLU_MISSING
#define TRIGGER_PATTERN_ADDRESS_0 EUDAQ_TLU_MISSING
#define TRIGGER_PATTERN_ADDRESS_1 EUDAQ_TLU_MISSING
#define AUX_PATTERN_ADDRESS_0 EUDAQ_TLU_MISSING
#define AUX_PATTERN_ADDRESS_1 EUDAQ_TLU_MISSING
#define ENABLE_DUT_VETO_ADDRESS EUDAQ_TLU_MISSING
#define TRIGGER_FSM_STATUS_VALUE_ADDRESS_0 EUDAQ_TLU_MISSING
#define TRIGGER_FSM_STATUS_VALUE_ADDRESS_1 EUDAQ_TLU_MISSING
#define TRIGGER_FSM_STATUS_VALUE_ADDRESS_2 EUDAQ_TLU_MISSING
#define CLOCK_GEN_RESET_BIT EUDAQ_TLU_MISSING
#define ENABLE_DMA_BIT EUDAQ_TLU_MISSING
#define RESET_DMA_COUNTER_BIT EUDAQ_TLU_MISSING

namespace tlu {

#define EUDAQ_TLU_REG(r) r,

  TLUAddresses v0_1 = {EUDAQ_TLU_REGISTERS 0};
}
