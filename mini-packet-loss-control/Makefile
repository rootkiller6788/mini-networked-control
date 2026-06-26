CC = gcc
CFLAGS = -std=c11 -Wall -Wextra -O2 -Iinclude
LDFLAGS = -lm
LIB = build/libpacket_loss_control.a
SRCS = $(wildcard src/*.c)
OBJS = $(patsubst src/%.c,build/%.o,$(SRCS))
TEST_SRC = tests/test_packet_loss.c
TEST_BIN = build/test_packet_loss
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
	$(CC) $(CFLAGS) -o $(TEST_BIN) $(TEST_SRC) -Lbuild -lpacket_loss_control $(LDFLAGS)
	./$(TEST_BIN)

examples: all $(EXAMPLE_BINS)
	@echo "All examples built."

build/%: examples/%.c $(LIB)
	$(CC) $(CFLAGS) -o $@ $< -Lbuild -lpacket_loss_control $(LDFLAGS)

demo: examples
	@echo "=== Running example1_bernoulli_channel ==="
	./build/example1_bernoulli_channel
	@echo ""
	@echo "=== Running example2_markov_channel ==="
	./build/example2_markov_channel
	@echo ""
	@echo "=== Running example3_inverted_pendulum ==="
	./build/example3_inverted_pendulum
	@echo ""
	@echo "=== Running example4_vehicle_platoon ==="
	./build/example4_vehicle_platoon

clean:
	rm -rf build
