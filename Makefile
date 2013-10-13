#Makefile

CC =gcc
LIB=-lm
FLAG=-g

SRC_DIR=.
BUILD_DIR=build
RESULT_DIR=.

NAME=main

RESULT=$(RESULT_DIR)/$(NAME)

HEADERS=\
		$(SRC_DIR)/util/util.h

OBJS=\
	 $(BUILD_DIR)/main.o\
	 $(BUILD_DIR)/bicgstab.o\
	 $(BUILD_DIR)/util.o

############################################################
#Build rules

all: dir $(RESULT)

dir:
	mkdir -p $(BUILD_DIR)

$(RESULT) : $(OBJS)
	$(CC) $(FLAG) $(OBJS) -o $(RESULT) $(LIB)

clean:
	rm -f $(OBJS) $(RESULT)
	rm -r $(BUILD_DIR)

############################################################
#individual files build  rules

$(BUILD_DIR)/main.o : $(SRC_DIR)/main.c $(HEADERS)
	$(CC) $(FLAG) -c $< -o $@

$(BUILD_DIR)/util.o : $(SRC_DIR)/util/util.c $(SRC_DIR)/util/util.h
	$(CC) $(FLAG) -c $< -o $@

$(BUILD_DIR)/bicgstab.o : $(SRC_DIR)/bicgstab/bicgstab.c $(SRC_DIR)/bicgstab/bicgstab.h
	$(CC) $(FLAG) -c $< -o $@

