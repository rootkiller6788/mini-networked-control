CC = gcc
CFLAGS = -std=c11 -Wall -Wextra -O2 -Iinclude
LDFLAGS = -lm
LIB = build/libtime_delay_system.a
SRCS = $(filter-out src/delay_applications.c,$(wildcard src/*.c))
OBJS = $(patsubst src/%.c,build/%.o,$(SRCS))
TEST_SRC = tests/test_time_delay.c
TEST_BIN = build/test_time_delay
EXAMPLE_SRCS = $(wildcard examples/*.c)
EXAMPLE_BINS = $(patsubst examples/%.c,build/%,$(EXAMPLE_SRCS))
APP_SRC = src/delay_applications.c
APP_BIN = build/delay_applications

.PHONY: all builddir clean test examples demo apps

all: builddir $(LIB)

builddir:
	@mkdir -p build

$(LIB): $(OBJS)
	$(AR) rcs $@ $^

build/%.o: src/%.c
	$(CC) $(CFLAGS) -c $< -o $@

test: all
	$(CC) $(CFLAGS) -o $(TEST_BIN) $(TEST_SRC) -Lbuild -ltime_delay_system $(LDFLAGS)
	./$(TEST_BIN)

examples: all $(EXAMPLE_BINS)
	@echo "All examples built."

build/%: examples/%.c $(LIB)
	$(CC) $(CFLAGS) -o $@ $< -Lbuild -ltime_delay_system $(LDFLAGS)

apps: all
	$(CC) $(CFLAGS) -o $(APP_BIN) $(APP_SRC) -Lbuild -ltime_delay_system $(LDFLAGS)
	./$(APP_BIN)

demo: examples apps
	@echo "=== Running example_smith_predictor ==="
	./build/example_smith_predictor
	@echo ""
	@echo "=== Running example_ncs_delay ==="
	./build/example_ncs_delay
	@echo ""
	@echo "=== Running example_dde_stability ==="
	./build/example_dde_stability

clean:
	rm -rf build
