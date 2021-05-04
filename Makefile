CFLAGS := -O3 -fPIC -g -std=c99 -pedantic -Wall -Wfatal-errors
INCLUDES := -Iinclude
LIBS := -Llib -lcube -Wl,-rpath,$(PWD)/lib -lm

.PHONY: lib clean

lib: lib/libcube.so
	@rm -f *.o

lib/libcube.so: src/cube.c include/cube.h
	@mkdir -p lib
	@$(CC) -o $@ $(CFLAGS) -shared $(INCLUDES) $<

clean:
	@rm -rf lib bin
