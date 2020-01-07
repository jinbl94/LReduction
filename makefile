TARGET=main
OBJS=main.o lalg.o
CC=g++
FLAGS= -lntl -g

default:$(TARGET)
$(TARGET):$(OBJS)
	$(CC) $(OBJS) -o $(TARGET) $(FLAGS)

%.o:%.cpp %.h
	$(CC) -c $< -o $@ $(FLAGS)

.PHONY:clean
clean:
	@rm -f $(OBJS) $(TARGET)
