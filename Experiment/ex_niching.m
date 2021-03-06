function ex_niching()
    test_algorithms = {
        {@DynamicNichingCSA, 1, 3, 1}, ...
        {@DynamicNichingCSA, 2, 3, 1}, ...
        {@DynamicNichingCSA, 3, 3, 1}, ...
        {@DynamicNichingCSA, 4, 3, 1}, ...
    };

    test_problems = generate_test_problems();
    
    running(test_algorithms, test_problems, "niching");
end