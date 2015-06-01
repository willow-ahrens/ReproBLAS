TARGETS := libindexedblas.a
SUBDIRS :=

INSTALL_LIB := $(TARGETS)

COGGED = damax.ccog damaxm.ccog                                        \
         zamax_sub.ccog zamaxm_sub.ccog                                \
         samax.ccog samaxm.ccog                                        \
         camax_sub.ccog camaxm_sub.ccog                                \
         dmdsum.ccog dmdasum.ccog dmdssq.ccog dmddot.ccog              \
         smssum.ccog smsasum.ccog smsssq.ccog smsdot.ccog              \
         cmcsum.ccog smcasum.ccog smcssq.ccog cmcdotu.ccog cmcdotc.ccog\
         zmzsum.ccog dmzasum.ccog dmzssq.ccog zmzdotu.ccog zmzdotc.ccog

PRECIOUS = damax.c damaxm.c                               \
           zamax_sub.c zamaxm_sub.c                       \
           samax.c samaxm.c                               \
           camax_sub.c camaxm_sub.c                       \
           smssum.c smsasum.c smsssq.c smsdot.c           \
           dmdsum.c dmdasum.c dmdssq.c dmddot.c           \
           cmcsum.c smcasum.c smcssq.c cmcdotu.c cmcdotc.c\
           zmzsum.c dmzasum.c dmzssq.c zmzdotu.c zmzdotc.c

LIBINDEXEDBLAS := $(OBJPATH)/libindexedblas.a

libindexedblas.a_DEPS = $$(LIBINDEXED)                                 \
                        damax.o damaxm.o                               \
                        zamax_sub.o zamaxm_sub.o                       \
                        samax.o samaxm.o                               \
                        camax_sub.o camaxm_sub.o                       \
                        dmdsum.o dmdasum.o dmdssq.o dmddot.o           \
                        smssum.o smsasum.o smsssq.o smsdot.o           \
                        cmcsum.o smcasum.o smcssq.o cmcdotu.o cmcdotc.o\
                        zmzsum.o dmzasum.o dmzssq.o zmzdotu.o zmzdotc.o\
                        didsum.o didasum.o didssq.o diddot.o           \
                        sissum.o sisasum.o sisssq.o sisdot.o           \
                        cicsum.o sicasum.o sicssq.o cicdotu.o cicdotc.o\
                        zizsum.o dizasum.o dizssq.o zizdotu.o zizdotc.o
                        #TODO dIAccum.o zIAccum.o sIAccum.o cIAccum.o        \
