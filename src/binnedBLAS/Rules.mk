TARGETS := libbinnedblas.a
SUBDIRS :=

INSTALL_LIB := $(TARGETS)

COGGED = damax.ccog damaxm.ccog                                         \
         zamax_sub.ccog zamaxm_sub.ccog                                 \
         samax.ccog samaxm.ccog                                         \
         camax_sub.ccog camaxm_sub.ccog                                 \
         dmdsum.ccog dmdasum.ccog dmdssq.ccog dmddot.ccog               \
         smssum.ccog smsasum.ccog smsssq.ccog smsdot.ccog               \
         cmcsum.ccog smcasum.ccog smcssq.ccog cmcdotu.ccog cmcdotc.ccog \
         zmzsum.ccog dmzasum.ccog dmzssq.ccog zmzdotu.ccog zmzdotc.ccog \
         dbdgemv.ccog dbdgemm.ccog                                      \
         sbsgemv.ccog sbsgemm.ccog                                      \
         cbcgemv.ccog cbcgemm.ccog                                      \
         zbzgemv.ccog zbzgemm.ccog                                      \

PRECIOUS = damax.c damaxm.c                                \
           zamax_sub.c zamaxm_sub.c                        \
           samax.c samaxm.c                                \
           camax_sub.c camaxm_sub.c                        \
           dmdsum.c dmdasum.c dmdssq.c dmddot.c            \
           smssum.c smsasum.c smsssq.c smsdot.c            \
           zmzsum.c dmzasum.c dmzssq.c zmzdotu.c zmzdotc.c \
           cmcsum.c smcasum.c smcssq.c cmcdotu.c cmcdotc.c \
           dbdgemv.c dbdgemm.c                             \
           sbsgemv.c sbsgemm.c                             \
           zbzgemv.c zbzgemm.c                             \
           cbcgemv.c cbcgemm.c                             \

LIBBINNEDBLAS := $(OBJPATH)/libbinnedblas.a

libbinnedblas.a_DEPS = $$(LIBBINNED)                                     \
                     damax.o damaxm.o                                \
                     zamax_sub.o zamaxm_sub.o                        \
                     samax.o samaxm.o                                \
                     camax_sub.o camaxm_sub.o                        \
                     dmdsum.o dmdasum.o dmdssq.o dmddot.o            \
                     zmzsum.o dmzasum.o dmzssq.o zmzdotu.o zmzdotc.o \
                     smssum.o smsasum.o smsssq.o smsdot.o            \
                     cmcsum.o smcasum.o smcssq.o cmcdotu.o cmcdotc.o \
                     dbdsum.o dbdasum.o dbdssq.o dbddot.o            \
                     zbzsum.o dbzasum.o dbzssq.o zbzdotu.o zbzdotc.o \
                     sbssum.o sbsasum.o sbsssq.o sbsdot.o            \
                     cbcsum.o sbcasum.o sbcssq.o cbcdotu.o cbcdotc.o \
                     dbdgemv.o dbdgemm.o                             \
                     zbzgemv.o zbzgemm.o                             \
                     sbsgemv.o sbsgemm.o                             \
                     cbcgemv.o cbcgemm.o

camax_sub.c_DEPS = camax_sub.ccog
camaxm_sub.c_DEPS = camaxm_sub.ccog
cbcgemm.c_DEPS = $$(GETTER) cbcgemm.ccog
cbcgemv.c_DEPS = $$(GETTER) cbcgemv.ccog
cmcdotc.c_DEPS = $$(GETTER) cmcdotc.ccog
cmcdotu.c_DEPS = $$(GETTER) cmcdotu.ccog
cmcsum.c_DEPS = $$(GETTER) cmcsum.ccog
damax.c_DEPS = damax.ccog
damaxm.c_DEPS = damaxm.ccog
dbdgemm.c_DEPS = $$(GETTER) dbdgemm.ccog
dbdgemv.c_DEPS = $$(GETTER) dbdgemv.ccog
dmdasum.c_DEPS = $$(GETTER) dmdasum.ccog
dmddot.c_DEPS = $$(GETTER) dmddot.ccog
dmdssq.c_DEPS = $$(GETTER) dmdssq.ccog
dmdsum.c_DEPS = $$(GETTER) dmdsum.ccog
dmzasum.c_DEPS = $$(GETTER) dmzasum.ccog
dmzssq.c_DEPS = $$(GETTER) dmzssq.ccog
samax.c_DEPS = samax.ccog
samaxm.c_DEPS = samaxm.ccog
sbsgemm.c_DEPS = $$(GETTER) sbsgemm.ccog
sbsgemv.c_DEPS = $$(GETTER) sbsgemv.ccog
smcasum.c_DEPS = $$(GETTER) smcasum.ccog
smcssq.c_DEPS = $$(GETTER) smcssq.ccog
smsasum.c_DEPS = $$(GETTER) smsasum.ccog
smsdot.c_DEPS = $$(GETTER) smsdot.ccog
smsssq.c_DEPS = $$(GETTER) smsssq.ccog
smssum.c_DEPS = $$(GETTER) smssum.ccog
zamax_sub.c_DEPS = zamax_sub.ccog
zamaxm_sub.c_DEPS = zamaxm_sub.ccog
zbzgemm.c_DEPS = $$(GETTER) zbzgemm.ccog
zbzgemv.c_DEPS = $$(GETTER) zbzgemv.ccog
zmzdotc.c_DEPS = $$(GETTER) zmzdotc.ccog
zmzdotu.c_DEPS = $$(GETTER) zmzdotu.ccog
zmzsum.c_DEPS = $$(GETTER) zmzsum.ccog
