CC = g++             
CFLAGS = -Wall -Werror -Wextra -O3 -std=c++17   
LIBRARY = s21_matrix_oop.a
SRC_FILES = s21_matrix_oop.cc
TEST_FILES = tests/s21_matrix_test.cc


OBJS = $(SRC_FILES:.cc=.o)
TEST_OBJ = $(TEST_FILES:.cc=.o)


all: $(LIBRARY) test
	rm $(OBJS)

test: CFLAGS += -lgtest -lgtest_main
test: $(TEST_OBJ) $(LIBRARY)
	$(CC) $(CFLAGS) -o s21_matrix_test $^
	./s21_matrix_test

$(LIBRARY): $(OBJS)
	ar rcs $@ $(OBJS)


%.o: %.cc
	$(CC) $(CFLAGS) -c $< -o $@

debug: 
	g++ main.cc s21_matrix_oop.cc 
	./a.out
	rm a.out

linters:
	cp ../materials/linters/.clang-format .
	clang-format -n *.cc *.h tests/*.cc
	rm .clang-format

clean:
	rm -f $(OBJS) $(LIBRARY)
	rm s21_matrix_test
	rm tests/s21_matrix_test.o