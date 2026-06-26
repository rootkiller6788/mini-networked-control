CC = gcc
CFLAGS = -std=c11 -Wall -Wextra -O2 -Iinclude
LDFLAGS = -lm
AR = ar

LIB = build/libcps_security.a
SRCS = $(wildcard src/*.c)
OBJS = $(patsubst src/%.c,build/%.o,$(SRCS))
TEST_SRC = tests/test_v2.c
TEST_BIN = build/test_cps_security
EXAMPLE_SRCS = $(wildcard examples/*.c)
EXAMPLE_BINS = $(patsubst examples/%.c,build/%,$(EXAMPLE_SRCS))
DEMO_SRC = demos/demo_v2.c
DEMO_BIN = build/demo_cps_security
BENCH_SRC = benches/bench_cps.c
BENCH_BIN = build/bench_cps

.PHONY: all builddir clean test examples demo bench

all: builddir $(LIB)

builddir:
	@mkdir -p build

$(LIB): $(OBJS)
	$(AR) rcs $@ $^

build/%.o: src/%.c
	$(CC) $(CFLAGS) -c $< -o $@

test: all
	$(CC) $(CFLAGS) -o $(TEST_BIN) $(TEST_SRC) -Lbuild -lcps_security $(LDFLAGS)
	./$(TEST_BIN)

examples: all $(EXAMPLE_BINS)
	@echo "All examples built."

build/%: examples/%.c $(LIB)
	$(CC) $(CFLAGS) -o $@ $< -Lbuild -lcps_security $(LDFLAGS)

demo: all
	$(CC) $(CFLAGS) -o $(DEMO_BIN) $(DEMO_SRC) -Lbuild -lcps_security $(LDFLAGS)
	./$(DEMO_BIN)

bench: all
	$(CC) $(CFLAGS) -o $(BENCH_BIN) $(BENCH_SRC) -Lbuild -lcps_security $(LDFLAGS)
	./$(BENCH_BIN)

run-examples: examples
	@echo "=== Running example_fdi_attack ==="
	./build/example_fdi_attack
	@echo ""
	@echo "=== Running example_replay_watermark ==="
	./build/example_replay_watermark
	@echo ""
	@echo "=== Running example_gametheory_defense ==="
	./build/example_gametheory_defense

clean:
	rm -rf build
