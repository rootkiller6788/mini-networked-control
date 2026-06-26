CC = gcc
CFLAGS = -std=c11 -Wall -Wextra -O2 -Iinclude
LDFLAGS = -lm
LIB = build/libbandwidth_limited_control.a
SRCS = $(wildcard src/*.c)
OBJS = $(patsubst src/%.c,build/%.o,$(SRCS))
TEST_SRC = tests/test_bandwidth_control.c
TEST_BIN = build/test_bandwidth_control
EXAMPLE_SRCS = $(wildcard examples/*.c)
EXAMPLE_BINS = $(patsubst examples/%.c,build/%,$(EXAMPLE_SRCS))

.PHONY: all builddir clean test examples demo

all: builddir $(LIB)

builddir:
	@mkdir -p build

$(LIB): $(OBJS)
	$(AR) rcs $@ $^

build/%.o: src/%.c
	$(CC) $(CFLAGS) -c $< -o $@

test: all
	$(CC) $(CFLAGS) -o $(TEST_BIN) $(TEST_SRC) -Lbuild -lbandwidth_limited_control $(LDFLAGS)
	./$(TEST_BIN)

examples: all $(EXAMPLE_BINS)
	@echo "All examples built."

build/%: examples/%.c $(LIB)
	$(CC) $(CFLAGS) -o $@ $< -Lbuild -lbandwidth_limited_control $(LDFLAGS)

demo: examples
	@echo "=== Running example1_bandwidth_limited_stabilization ==="
	./build/example1_bandwidth_limited_stabilization
	@echo ""
	@echo "=== Running example2_event_triggered_control ==="
	./build/example2_event_triggered_control
	@echo ""
	@echo "=== Running example3_multi_loop_scheduling ==="
	./build/example3_multi_loop_scheduling

clean:
	rm -rf build
