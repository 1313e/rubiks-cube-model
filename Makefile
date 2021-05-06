CFLAGS := -O3 -fPIC -g -std=c99 -pedantic -Wall -Wfatal-errors
INCLUDES := -Iinclude
LIBS := -Llib -lm

.PHONY: lib clean

lib: lib/libcube.so
	@rm -f *.o

lib/libcube.so: src/cube.c include/cube.h
	@mkdir -p lib
	@$(CC) -o $@ $(CFLAGS) -shared $(INCLUDES) $(LIBS) $<

clean:
	@rm -rf lib bin
