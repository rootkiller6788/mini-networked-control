CC = gcc
CFLAGS = -std=c11 -Wall -Wextra -O2 -Iinclude
LDFLAGS = -lm
LIB = build/libebc.a
SRCS = $(wildcard src/*.c)
OBJS = $(patsubst src/%.c,build/%.o,$(SRCS))
TEST_SRCS = $(wildcard tests/*.c)
TEST_BINS = $(patsubst tests/%.c,build/test_%,$(TEST_SRCS))
EXAMPLE_SRCS = $(wildcard examples/*.c)
EXAMPLE_BINS = $(patsubst examples/%.c,build/example_%,$(EXAMPLE_SRCS))
all: builddir $(LIB)
builddir: ; @mkdir -p build
$(LIB): $(OBJS) ; $(AR) rcs $@ $^
build/%.o: src/%.c | builddir ; $(CC) $(CFLAGS) -c $< -o $@
build/test_%: tests/%.c $(LIB) ; $(CC) $(CFLAGS) $< -Lbuild -lebc $(LDFLAGS) -o $@
build/example_%: examples/%.c $(LIB) ; $(CC) $(CFLAGS) $< -Lbuild -lebc $(LDFLAGS) -o $@
test: $(LIB) $(TEST_BINS) ; for t in $(TEST_BINS); do $$t || exit 1; done && echo All tests passed
examples: $(LIB) $(EXAMPLE_BINS) ; for e in $(EXAMPLE_BINS); do $$e; done
clean: ; rm -rf build
