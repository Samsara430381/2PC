macro (add_test_STR _name)
	add_executable(${_name}-STR "test_${_name}.cpp")
    target_link_libraries(${_name}-STR SCI-OT) 
endmacro()


add_test_STR(Exp_N)
# add_test_STR(Exp_P)
add_test_STR(exp)
# add_test_STR(NewtonRaphson)
add_test_STR(Sigmoid)
add_test_STR(Tanh)