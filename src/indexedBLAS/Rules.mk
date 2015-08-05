TARGETS := libindexedblas.a
SUBDIRS :=

INSTALL_LIB := $(TARGETS)

COGGED = damax.ccog damaxm.ccog                                        \
         zamax_sub.ccog zamaxm_sub.ccog                                \
         samax.ccog samaxm.ccog                                        \
         camax_sub.ccog camaxm_sub.ccog                                \
         dmdsum.ccog dmdasum.ccog dmdssq.ccog dmddot.ccog dmdgemv.ccog \
         smssum.ccog smsasum.ccog smsssq.ccog smsdot.ccog              \
         cmcsum.ccog smcasum.ccog smcssq.ccog cmcdotu.ccog cmcdotc.ccog\
         zmzsum.ccog dmzasum.ccog dmzssq.ccog zmzdotu.ccog zmzdotc.ccog zmzgemv.ccog

PRECIOUS = damax.c damaxm.c                               \
           zamax_sub.c zamaxm_sub.c                       \
           samax.c samaxm.c                               \
           camax_sub.c camaxm_sub.c                       \
           dmdsum.c dmdasum.c dmdssq.c dmddot.c dmdgemv.c \
           smssum.c smsasum.c smsssq.c smsdot.c           \
           zmzsum.c dmzasum.c dmzssq.c zmzdotu.c zmzdotc.c zmzgemv.c\
           cmcsum.c smcasum.c smcssq.c cmcdotu.c cmcdotc.c\

LIBINDEXEDBLAS := $(OBJPATH)/libindexedblas.a

libindexedblas.a_DEPS = $$(LIBINDEXED)                                 \
                        damax.o damaxm.o                               \
                        zamax_sub.o zamaxm_sub.o                       \
                        samax.o samaxm.o                               \
                        camax_sub.o camaxm_sub.o                       \
                        dmdsum.o dmdasum.o dmdssq.o dmddot.o dmdgemv.o \
                        smssum.o smsasum.o smsssq.o smsdot.o           \
                        cmcsum.o smcasum.o smcssq.o cmcdotu.o cmcdotc.o\
                        zmzsum.o dmzasum.o dmzssq.o zmzdotu.o zmzdotc.o zmzgemv.o\
                        didsum.o didasum.o didssq.o diddot.o didgemv.o \
                        sissum.o sisasum.o sisssq.o sisdot.o           \
                        cicsum.o sicasum.o sicssq.o cicdotu.o cicdotc.o\
                        zizsum.o dizasum.o dizssq.o zizdotu.o zizdotc.o zizgemv.o\

damax.c_DEPS = damax.ccog
damaxm.c_DEPS = damaxm.ccog
zamax_sub.c_DEPS = zamax_sub.ccog
zamaxm_sub.c_DEPS = zamaxm_sub.ccog
samax.c_DEPS = samax.ccog
samaxm.c_DEPS = samaxm.ccog
camax_sub.c_DEPS = camax_sub.ccog
camaxm_sub.c_DEPS = camaxm_sub.ccog
dmdsum.c_DEPS = $$(GETTER) dmdsum.ccog
dmdasum.c_DEPS = $$(GETTER) dmdasum.ccog
dmdssq.c_DEPS = $$(GETTER) dmdssq.ccog
dmddot.c_DEPS = $$(GETTER) dmddot.ccog
dmdgemv.c_DEPS = $$(GETTER) dmdgemv.ccog
smssum.c_DEPS = $$(GETTER) smssum.ccog
smsasum.c_DEPS = $$(GETTER) smsasum.ccog
smsssq.c_DEPS = $$(GETTER) smsssq.ccog
smsdot.c_DEPS = $$(GETTER) smsdot.ccog
cmcsum.c_DEPS = $$(GETTER) cmcsum.ccog
smcasum.c_DEPS = $$(GETTER) smcasum.ccog
smcssq.c_DEPS = $$(GETTER) smcssq.ccog
cmcdotu.c_DEPS = $$(GETTER) cmcdotu.ccog
cmcdotc.c_DEPS = $$(GETTER) cmcdotc.ccog
zmzsum.c_DEPS = $$(GETTER) zmzsum.ccog
dmzasum.c_DEPS = $$(GETTER) dmzasum.ccog
dmzssq.c_DEPS = $$(GETTER) dmzssq.ccog
zmzdotu.c_DEPS = $$(GETTER) zmzdotu.ccog
zmzdotc.c_DEPS = $$(GETTER) zmzdotc.ccog
zmzgemv.c_DEPS = $$(GETTER) zmzgemv.ccog
