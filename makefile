CC=gcc -g
PP=/home/glenn/lib/penplot.a
FC=g77 -O

ch2: ch2.o cosft.o sinft.o realftshift.o four1shift.o realft.o four1.o
	${CC} -o $@ ch2.o cosft.o sinft.o realftshift.o four1shift.o realft.o four1.o -lm
ch2.o:ch2.c ch2.h init.c

ch2mu: ch2mu.o cosft.o sinft.o realftshift.o four1shift.o realft.o four1.o
	${CC} -o $@ ch2mu.o cosft.o sinft.o realftshift.o four1shift.o realft.o four1.o -lm
ch2mu.o:ch2mu.c ch2.h init.c

cjp2: cjp2.o cosft.o sinft.o realftshift.o four1shift.o realft.o four1.o
	${CC} -o $@ cjp2.o cosft.o sinft.o realftshift.o four1shift.o realft.o four1.o -lm
channelphys: channelphys.o cosft.o sinft.o realftshift.o four1shift.o realft.o four1.o
	${CC} -o $@ channelphys.o cosft.o sinft.o realftshift.o four1shift.o realft.o four1.o -lm
cjp2.o:cjp2.c
channelphys.o:channelphys.c
cosft.o:cosft.c
sinft.o:sinft.c
realftshift.o:realftshift.c
four1shift.o:four1shift.c
realft.o:realft.c
four1.o:four1.c

cj2:cj2.c
	gcc -O -o cj2 cj2.c sinft.o cosft.o realftshift.o four1shift.o realft.o four1.o -lm
cj2r:cj2r.c
	gcc -O -o $@ $@.c sinft.o cosft.o realftshift.o four1shift.o realft.o four1.o -lm
conv: conv.o phys.o transforms.o ../cosft.o ../sinft.o ../realftshift.o ../four1shift.o ../realft.o ../four1.o
	${CC} conv.o phys.o transforms.o ../cosft.o ../sinft.o ../realftshift.o ../four1shift.o ../realft.o ../four1.o -O -lm -o $@

convm: convm.o transforms.o ../cosft.o ../sinft.o ../realftshift.o ../four1shift.o ../realft.o ../four1.o
	${CC} convm.o transforms.o ../cosft.o ../sinft.o ../realftshift.o ../four1shift.o ../realft.o ../four1.o -O -lm -o $@
plotitm:plotitm.f ../icread.o
	${FC} plotitm.f ../icread.o ${PP} -o $@

ctest: ctest.o ../transforms.o ../cosft.o ../sinft.o ../realftshift.o ../four1shift.o ../realft.o ../four1.o
	${CC} ctest.o ../transforms.o ../cosft.o ../sinft.o ../realftshift.o ../four1shift.o ../realft.o ../four1.o -O -lm -o $@

plotit:plotit.f ../icread.o
	${FC} plotit.f ../icread.o ${PP} -o $@

tracer: tracer.o trphys.o transforms.o cosft.o sinft.o realftshift.o four1shift.o realft.o four1.o
	${CC} tracer.o trphys.o transforms.o cosft.o sinft.o realftshift.o four1shift.o realft.o four1.o -O -lm -o $@

trplotit:trplotit.f icread.o
	${FC} trplotit.f icread.o ${PP} -o $@

trace3: trace3.f
	${FC} trace3.f /home/glenn/lib/penplot.a -o $@
