TARGETS := libindexedblas.a
SUBDIRS := 

INSTALL_LIB := $(TARGETS)

COGGED = dasumI2.ccog dsumI2.ccog dnrm2I2.ccog ddotI2.ccog                 \
         zsumI2.ccog dzasumI2.ccog dznrm2I2.ccog zdotcI2.ccog zdotuI2.ccog \
         sasumI2.ccog ssumI2.ccog snrm2I2.ccog sdotI2.ccog                 \
         csumI2.ccog scasumI2.ccog scnrm2I2.ccog cdotcI2.ccog cdotuI2.ccog \

PRECIOUS = dasumI2.c dsumI2.c dnrm2I2.c ddotI2.c              \
           zsumI2.c dzasumI2.c dznrm2I2.c zdotcI2.c zdotuI2.c \
           sasumI2.c ssumI2.c snrm2I2.c sdotI2.c              \
           csumI2.c scasumI2.c scnrm2I2.c cdotcI2.c cdotuI2.c \

LIBINDEXEDBLAS := $(OBJPATH)/libindexedblas.a

libindexedblas.a_DEPS = $$(LIBINDEXED)                                     \
                        dasumI2.o dsumI2.o dnrm2I2.o ddotI2.o              \
                        zsumI2.o dzasumI2.o dznrm2I2.o zdotcI2.o zdotuI2.o \
                        sasumI2.o ssumI2.o snrm2I2.o sdotI2.o              \
                        csumI2.o scasumI2.o scnrm2I2.o cdotcI2.o cdotuI2.o \
                        dsumI.o dasumI.o ddotI.o dnrm2I.o                  \
                        zdotI.o zsumI.o dzasumI.o dznrm2I.o                \
                        ssumI.o sasumI.o sdotI.o snrm2I.o                  \
                        csumI.o scasumI.o cdotI.o scnrm2I.o                \
                        damax.o damaxm.o                                   \
                        zamax.o zamaxm.o                                   \
                        samax.o samaxm.o                                   \
                        camax.o camaxm.o                                   \
                        dIAccum.o zIAccum.o sIAccum.o cIAccum.o
