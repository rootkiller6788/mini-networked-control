CC = gcc
CFLAGS = -std=c11 -Wall -Wextra -O2 -Iinclude
LDFLAGS = -lm
LIB = build/libconsensus.a
SRCS = $(wildcard src/*.c)
OBJS = $(patsubst src/%.c,build/%.o,$(SRCS))
TEST_SRC = tests/test_consensus.c
TEST_BIN = build/test_consensus
EXAMPLE_SRCS = $(wildcard examples/*.c)
EXAMPLE_BINS = $(patsubst examples/%.c,build/%,$(EXAMPLE_SRCS))

.PHONY: all builddir clean test examples demo

all: builddir $(LIB)

builddir:
	mkdir -p build

$(LIB): $(OBJS)
	$(AR) rcs $@ $^

build/%.o: src/%.c
	$(CC) $(CFLAGS) -c $< -o $@

test: all
	$(CC) $(CFLAGS) -o $(TEST_BIN) $(TEST_SRC) -Lbuild -lconsensus $(LDFLAGS)
	./$(TEST_BIN)

examples: all $(EXAMPLE_BINS)
	@echo "All examples built."

build/%: examples/%.c $(LIB)
	$(CC) $(CFLAGS) -o $@ $< -Lbuild -lconsensus $(LDFLAGS)

demo: examples
	@echo "=== Running example_averaging ==="
	./build/example_averaging
	@echo ""
	@echo "=== Running example_rendezvous ==="
	./build/example_rendezvous
	@echo ""
	@echo "=== Running example_clock_sync ==="
	./build/example_clock_sync
	@echo ""
	@echo "=== Running example_formation ==="
	./build/example_formation

clean:
	rm -rf build