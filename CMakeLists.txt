cmake_minimum_required(VERSION 3.9.6)
project(notes_extracting)

set(CMAKE_CXX_STANDARD 11)

include_directories(./)
include_directories(./MFE)
include_directories(./ScoreFollowing)


set(TOP_FILES wave.c main.cpp)
set(NMF_FILES MFE/MFE.c MFE/resample.cpp MFE/kiss_fft.c)
file(GLOB SCORE_FILES ScoreFollowing/*.cpp)
list(REMOVE_ITEM SCORE_FILES ${CMAKE_CURRENT_SOURCE_DIR}/ScoreFollowing/main.cpp)

add_executable(notes_extracting ${NMF_FILES} ${SCORE_FILES} ${TOP_FILES})