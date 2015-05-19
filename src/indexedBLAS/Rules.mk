TARGETS := libindexedblas.a
SUBDIRS :=

INSTALL_LIB := $(TARGETS)

COGGED = damax.ccog damaxm.ccog                                        \
         zamax_sub.ccog zamaxm_sub.ccog                                \
         samax.ccog samaxm.ccog                                        \
         camax_sub.ccog camaxm_sub.ccog                                \
         dmdsum.ccog dmdasum.ccog dmdnrm.ccog dmddot.ccog              \
         smssum.ccog smsasum.ccog smsnrm.ccog smsdot.ccog              \
         cmcsum.ccog smcasum.ccog smcnrm.ccog cmcdotu.ccog cmcdotc.ccog\
         zmzsum.ccog dmzasum.ccog dmznrm.ccog zmzdotu.ccog zmzdotc.ccog

PRECIOUS = damax.c damaxm.c                               \
           zamax_sub.c zamaxm_sub.c                       \
           samax.c samaxm.c                               \
           camax_sub.c camaxm_sub.c                       \
           smssum.c smsasum.c smsnrm.c smsdot.c           \
           dmdsum.c dmdasum.c dmdnrm.c dmddot.c           \
           cmcsum.c smcasum.c smcnrm.c cmcdotu.c cmcdotc.c\
           zmzsum.c dmzasum.c dmznrm.c zmzdotu.c zmzdotc.c

LIBINDEXEDBLAS := $(OBJPATH)/libindexedblas.a

libindexedblas.a_DEPS = $$(LIBINDEXED)                                 \
                        damax.o damaxm.o                               \
                        zamax_sub.o zamaxm_sub.o                       \
                        samax.o samaxm.o                               \
                        camax_sub.o camaxm_sub.o                       \
                        dmdsum.o dmdasum.o dmdnrm.o dmddot.o           \
                        smssum.o smsasum.o smsnrm.o smsdot.o           \
                        cmcsum.o smcasum.o smcnrm.o cmcdotu.o cmcdotc.o\
                        zmzsum.o dmzasum.o dmznrm.o zmzdotu.o zmzdotc.o\
                        didsum.o didasum.o didnrm.o diddot.o           \
                        sissum.o sisasum.o sisnrm.o sisdot.o           \
                        cicsum.o sicasum.o sicnrm.o cicdotu.o cicdotc.o\
                        zizsum.o dizasum.o diznrm.o zizdotu.o zizdotc.o
#                        dgemvI.o
                        #TODO rename *_sub files
                        #TODO dIAccum.o zIAccum.o sIAccum.o cIAccum.o        \
