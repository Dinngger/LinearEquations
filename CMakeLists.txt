project(LinearEquations)

cmake_minimum_required(VERSION 3.10)

aux_source_directory(. DIR_SRCS)

add_executable(LinearEquationsTest ${DIR_SRCS})
