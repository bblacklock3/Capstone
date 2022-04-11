prop = load_prop('prototype');
prop.max_hip_torque = 0.8;
prop.max_knee_torque = 0.5;
tic
infant_id([1 1],[1 1],[1 1],prop)
toc