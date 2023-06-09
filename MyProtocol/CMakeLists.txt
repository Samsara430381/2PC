add_library(SCI-MyProtocol MyProtocol.cpp)
target_link_libraries(SCI-MyProtocol
    PUBLIC SCI-LinearOT SCI-BuildingBlocks
)
