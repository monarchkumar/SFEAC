all: makebin main #version

makebin:
	@mkdir -p bin

main:
	@g++ ./src/main.cc -o ./bin/main -std=c++20

# version:
# 	g++ ./src/version.cc -o ./bin/version -std=c++20

clean:
	@rm ./bin/*
