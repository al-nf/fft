CC := g++
CFLAGS := -Wall -Wextra -pedantic

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
    CC := g++-14
endif

all:
	$(CC) $(CFLAGS) main.cpp -o fft

clean:
	rm -f fft
