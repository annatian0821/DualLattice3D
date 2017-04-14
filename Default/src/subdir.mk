################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/ConnectLat.cpp \
../src/Dataclass.cpp \
../src/FluidCal.cpp \
../src/FluidNetwork.cpp \
../src/GConjGrad.cpp \
../src/GLatRemoval.cpp \
../src/InitCond.cpp \
../src/Input.cpp \
../src/Loading.cpp \
../src/MiscTools.cpp \
../src/Output.cpp \
../src/PCG_Eigen.cpp \
../src/ProProcess.cpp \
../src/SparseMatrix.cpp \
../src/Statistics.cpp \
../src/Stress.cpp \
../src/Tesellation.cpp \
../src/Testing.cpp \
../src/Vertices.cpp \
../src/main.cpp 

OBJS += \
./src/ConnectLat.o \
./src/Dataclass.o \
./src/FluidCal.o \
./src/FluidNetwork.o \
./src/GConjGrad.o \
./src/GLatRemoval.o \
./src/InitCond.o \
./src/Input.o \
./src/Loading.o \
./src/MiscTools.o \
./src/Output.o \
./src/PCG_Eigen.o \
./src/ProProcess.o \
./src/SparseMatrix.o \
./src/Statistics.o \
./src/Stress.o \
./src/Tesellation.o \
./src/Testing.o \
./src/Vertices.o \
./src/main.o 

CPP_DEPS += \
./src/ConnectLat.d \
./src/Dataclass.d \
./src/FluidCal.d \
./src/FluidNetwork.d \
./src/GConjGrad.d \
./src/GLatRemoval.d \
./src/InitCond.d \
./src/Input.d \
./src/Loading.d \
./src/MiscTools.d \
./src/Output.d \
./src/PCG_Eigen.d \
./src/ProProcess.d \
./src/SparseMatrix.d \
./src/Statistics.d \
./src/Stress.d \
./src/Tesellation.d \
./src/Testing.d \
./src/Vertices.d \
./src/main.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O2 -g -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


