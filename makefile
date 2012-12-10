VPATH = ./src
CC = g++
CFLAGS = -Wall -O3


all : engine parser plot

clear:
	rm engine
	rm parser
	rm plot

remake:
	rm engine
	rm parser
	make parser
	./parser
	make photo

engine : engine.cpp
	$(CC) $(CFLAGS) -o $@ $^

plot : plot.cpp
	$(CC) $(CFLAGS) -o $@ -lpng -lgd $^

parser : parser.cpp
	$(CC) $(CFLAGS) -o $@ $^

