function ex_dynamic()
    test_algorithms = {
        {@DynamicNichingCSA, 3, 1, 1}, ...
        {@DynamicNichingCSA, 3, 2, 1}, ...
        {@DynamicNichingCSA, 3, 3, 1}, ...
    };

    test_problems = generate_test_problems();
    
    running(test_algorithms, test_problems, "dynamic");
end