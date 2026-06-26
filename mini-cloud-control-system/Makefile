CC = gcc
CFLAGS = -std=c11 -Wall -Wextra -O2 -Iinclude
LDFLAGS = -lm
AR = ar
LIB = build/libcloud_control.a
SRCS = $(wildcard src/*.c)
OBJS = $(patsubst src/%.c,build/%.o,$(SRCS))
TEST_SRC = tests/test_cloud_control.c
TEST_BIN = build/test_cloud_control
EXAMPLE_SRCS = $(wildcard examples/*.c)
EXAMPLE_BINS = $(patsubst examples/%.c,build/%,$(EXAMPLE_SRCS))
BENCH_SRC = benches/bench_cloud_control.c
BENCH_BIN = build/bench_cloud_control
DEMO_SRC = demos/demo_overview.c
DEMO_BIN = build/demo_overview

.PHONY: all builddir clean test examples demo bench

all: builddir $(LIB)

builddir:
	@mkdir -p build

$(LIB): $(OBJS)
	$(AR) rcs $@ $^

build/%.o: src/%.c
	$(CC) $(CFLAGS) -c $< -o $@

test: all
	$(CC) $(CFLAGS) -o $(TEST_BIN) $(TEST_SRC) -Lbuild -lcloud_control $(LDFLAGS)
	./$(TEST_BIN)

examples: all $(EXAMPLE_BINS)
	@echo "All examples built."

build/%: examples/%.c $(LIB)
	$(CC) $(CFLAGS) -o $@ $< -Lbuild -lcloud_control $(LDFLAGS)

demo: examples
	@echo "=== Running Cloud Control Demo ==="
	./build/demo_cloud_control_basic
	@echo ""
	@echo "=== Running Delay Compensation Example ==="
	./build/demo_delay_compensation
	@echo ""
	@echo "=== Running Resource Scaling Example ==="
	./build/demo_resource_scaling

bench: all
	$(CC) $(CFLAGS) -o $(BENCH_BIN) $(BENCH_SRC) -Lbuild -lcloud_control $(LDFLAGS)
	./$(BENCH_BIN)

clean:
	rm -rf build
