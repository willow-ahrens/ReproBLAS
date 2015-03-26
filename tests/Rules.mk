TARGETS := libindexedblas.a
SUBDIRS := 

INSTALL_LIB := $(TARGETS)

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
