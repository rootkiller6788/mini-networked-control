CC = gcc
CFLAGS = -Wall -Wextra -O2 -std=c99 -Iinclude
LDFLAGS = -lm

SRC = $(wildcard src/*.c)
OBJ = $(SRC:.c=.o)
TARGET = libquantizedcontrol.a
TEST_SRC = tests/test_qc.c
TEST_BIN = test_qc

.PHONY: all clean test check

all: $(TARGET)

$(TARGET): $(OBJ)
	ar rcs $@ $^

src/%.o: src/%.c
	$(CC) $(CFLAGS) -c $< -o $@

test: $(TARGET) $(TEST_SRC)
	$(CC) $(CFLAGS) $(TEST_SRC) -o $(TEST_BIN) -L. -lquantizedcontrol $(LDFLAGS)
	./$(TEST_BIN)

clean:
	rm -f src/*.o $(TARGET) $(TEST_BIN) $(TEST_BIN).exe

check: test
	@echo "All tests passed."
