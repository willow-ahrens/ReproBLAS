#ifndef _UNIT_WRAPPER_H
#define _UNIT_WRAPPER_H

#define UNIT_FLOPS 0
#define UNIT_HERTZ 1

const char* unit_name(int unit){
  switch(unit){
    case UNIT_FLOPS:
      return "FLOPS";
    case UNIT_HERTZ:
      return "Hertz";
  }
  return "";
}

#endif
