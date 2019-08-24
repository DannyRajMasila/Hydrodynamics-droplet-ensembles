% Objective function as described in the manuscript is written out for
% selected shapes. For custom shapes, the user may have to input the
% necessary OBJfun based on the idea used here.

switch shape
    case 'circle5'
        ObjFn =[0.00022
            0.0002363
            0.0002539
            0.0002135
            0.0002355];
    case 'circle6'
        ObjFn =1.0e-03*[
            0.231811791237006
            0.236575987343484
            0.252427238148078
            0.252575505595799
            0.253631562809180
            0.266515669295576];
    case 'square8'
        ObjFn=[2.2138087e-04
            2.2660178e-04
            2.2751093e-04
            2.3362750e-04
            3.1013959e-04
            3.1461595e-04
            3.1606743e-04
            3.1730420e-04];
    case 'triangle6'
        ObjFn =1.0e-03*[0.1503
            0.1530
            0.1547
            0.3016
            0.3077
            0.3239];
end