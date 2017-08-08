function [X, Z] = eskf_map(sensor, varargin)
% ESKF Error State Kalman Filter with IMU as inputs
%
% INPUTS:
%   sensor - struct stored in provided dataset, fields include
%          - is_ready: logical, indicates whether sensor data is valid
%          - t: sensor timestamp
%          - rpy, omg, acc: imu readings
%          - img: uint8, 240x376 grayscale image
%          - id: 1xn ids of detected tags
%          - p0, p1, p2, p3, p4: 2xn pixel position of center and
%                                four corners of detected tags
%            Y
%            ^ P3 == P2
%            | || P0 ||
%            | P4 == P1
%            o---------> X
%   varargin - any variables you wish to pass into the function, could be
%              a data structure to represent the map or camera parameters,
%              your decision. But for the purpose of testing, since we don't
%              know what inputs you will use, you have to specify them in
%              init_script by doing
%              ekf1_handle = ...
%                  @(sensor) ekf2(sensor, your input arguments);
%
% OUTPUTS:
% X - nx1 state of the quadrotor, n should be greater or equal to 9
%     the state should be in the following order
%     [x; y; z; vx; vy; vz; qw; qx; qy; qz; other states you use]
%     we will only take the first 10 rows of X
%% Initialization
R_b_c = [sqrt(2) / 2, -sqrt(2) / 2, 0;
        -sqrt(2) / 2, -sqrt(2) / 2, 0;
         0          , 0           ,-1];
R_c_b = R_b_c';
T_b_c = [-0.04;0;-0.03];
T_c_b = -R_c_b * T_b_c;
K = [314.1779 0         199.4848;
     0        314.2218  113.7838;
     0        0         1       ];
% We are using the second method here.
% What we want to do here is to separate the robot state and the features.
% As in the phase 3, the eskf2 works very well.
persistent sensor_history;
persistent p_history;
persistent q_history;
persistent v_history;
persistent bias_acceleration_history;
persistent bias_gyro_history;
persistent mu;
persistent sigma;
persistent feature_history;
persistent id_history;
persistent feature_sigma;
if isempty(sensor_history)
    dt = sensor.t;
    p_history = [0;0;0];
    q_history = [1;0;0;0];
    v_history = [0;0;0];
    bias_acceleration_history = [0;0;0];
    bias_gyro_history = [0;0;0];
    mu = zeros(15,1);
    sigma = 0.01 * eye(15);
    sensor_history = sensor;
    feature_history = [];
    feature_sigma = [];
    
else
    dt = sensor.t - sensor_history.t;
end
g = [0;0;-9.8];
R = quat2rot(q_history);
omega_imu = sensor.omg;
a_imu = sensor.acc;
% Q = [0.01 * eye(3),     zeros(3,3),         zeros(3,3),         zeros(3,3);
%      zeros(3,3),        eye(3) * dt^2 * 10, zeros(3,3),         zeros(3,3);
%      zeros(3,3),        zeros(3,3),         0.01 * eye(3),      zeros(3,3);
%      zeros(3,3),        zeros(3,3),         zeros(3,3),         eye(3) * dt * 0.1];
Q = zeros(12,12);

%% Prediction
v_prediction = v_history + dt * (R * (a_imu - bias_acceleration_history) + g);
p_prediction = p_history + v_history * dt + 0.5 * dt^2 * (R * (a_imu - bias_acceleration_history) + g);
d_theta = (omega_imu - bias_gyro_history) * dt;
k = d_theta / norm(d_theta);
theta = norm(d_theta);
dq_prediction = [cos(theta / 2); k(1) * sin(theta / 2); k(2) * sin(theta / 2); k(3) * sin(theta / 2)];
q_prediction = quatmultiply(q_history, dq_prediction);
R_q_prediction = quat2rot(dq_prediction);
bias_acceleration_prediction = bias_acceleration_history;
bias_gyro_prediction = bias_gyro_history;
Fx = [eye(3),     eye(3) * dt, zeros(3,3),                                        zeros(3,3), zeros(3,3);
      zeros(3,3), eye(3),      -R * skew(a_imu - bias_acceleration_history) * dt, -R * dt,    zeros(3,3);
      zeros(3,3), zeros(3,3),  R_q_prediction',                                   zeros(3,3)  -eye(3) * dt;
      zeros(3,3), zeros(3,3),  zeros(3,3),                                        eye(3),     zeros(3,3)
      zeros(3,3), zeros(3,3),  zeros(3,3),                                        zeros(3,3)  eye(3)];
Fi = [zeros(3,3), zeros(3,3), zeros(3,3), zeros(3,3);
      eye(3),     zeros(3,3), zeros(3,3), zeros(3,3);
      zeros(3,3), eye(3),     zeros(3,3), zeros(3,3);
      zeros(3,3), zeros(3,3), eye(3),     zeros(3,3);
      zeros(3,3), zeros(3,3), zeros(3,3), eye(3)];
mu_bar = Fx * mu;
sigma_bar = Fx * sigma * Fx' + Fi * Q * Fi';
T_b_w = p_history;
R_b_w = quat2rot(q_history);
sensor_history = sensor;
%% Feature Sensing
if ~isempty(sensor.id)
    new_tag_id = sensor.id(~ismember(sensor.id, id_history));
    id_history = [id_history, new_tag_id];
    if ~isempty(new_tag_id)
        new_tag_number = size(new_tag_id, 2);
        new_tag_world_coordinate = zeros(3, 4 * new_tag_number);
        for i = 1:new_tag_number
            new_tag_world_coordinate(1:3, 4 * i - 3) = ...
            R_b_w * (R_c_b * K \ [sensor.p1(:, sensor.id == new_tag_id(i)); 1] + T_c_b) + T_b_w;
            new_tag_world_coordinate(1:3, 4 * i - 2) = ...
            R_b_w * (R_c_b * K \ [sensor.p2(:, sensor.id == new_tag_id(i)); 1] + T_c_b) + T_b_w;
            new_tag_world_coordinate(1:3, 4 * i - 1) = ...
            R_b_w * (R_c_b * K \ [sensor.p3(:, sensor.id == new_tag_id(i)); 1] + T_c_b) + T_b_w;
            new_tag_world_coordinate(1:3, 4 * i) = ...
            R_b_w * (R_c_b * K \ [sensor.p4(:, sensor.id == new_tag_id(i)); 1] + T_c_b) + T_b_w;
        end
        feature_expansion = reshape(new_tag_world_coordinate, 12 * new_tag_number, 1);
        feature_history = [feature_history; feature_expansion];
        feature_sigma = blkdiag(feature_sigma, 1000 * eye(12 * new_tag_number));
    end
end

%% Update
if isempty(sensor.id)
    p_history = p_prediction;
    q_history = q_prediction;
    v_history = v_prediction;
    bias_acceleration_history = bias_acceleration_prediction;
    bias_gyro_history = bias_gyro_prediction;
    mu = mu_bar;
    sigma = sigma_bar;
else
    sigma_together = blkdiag(sigma_bar, feature_sigma);
    tag_number_seen = size(sensor.id,2);
    actual_measurement = [zeros(2, 4 * tag_number_seen); ones(1, 4 * tag_number_seen)];
    camera_coordinate = zeros(3, 4 * tag_number_seen);
    predicted_measurement = zeros(2, 4 * tag_number_seen);
    C = zeros(8 * tag_number_seen, 15 + size(feature_history, 1));
    for j = 1 : tag_number_seen
        actual_measurement(1:2, 4 * j - 3) = sensor.p1(:, j);
        actual_measurement(1:2, 4 * j - 2) = sensor.p2(:, j);
        actual_measurement(1:2, 4 * j - 1) = sensor.p3(:, j);
        actual_measurement(1:2, 4 * j) = sensor.p4(:, j);
        tag_location = find(id_history == sensor.id(j), 1);
        p1_world_coordinate = feature_history(12 * (tag_location - 1) + 1 : 12 * (tag_location - 1) + 3, 1);
        camera_coordinate(:, 4 * j - 3) = R_b_c * (R_b_w' * p1_world_coordinate - R_b_w' * T_b_w) - T_b_c;
        X1 = camera_coordinate(1, 4 * j - 3);
        Y1 = camera_coordinate(2, 4 * j - 3);
        Z1 = camera_coordinate(3, 4 * j - 3);
        predicted_measurement(:, 4 * j - 3) = [X1/Z1; Y1/Z1];
        Jp_1 = [1, 0, -X1/Z1; 0, 1, -Y1/Z1] / Z1;
        C(8 * (j - 1) + 1 : 8 * (j - 1) + 2, 1 : 3) = -Jp_1 * R_b_c * R_b_w';
        C(8 * (j - 1) + 1 : 8 * (j - 1) + 2, 7 : 9) = Jp_1 * R_b_c * R_b_w' * skew(p1_world_coordinate - T_b_w);
        C(8 * (j - 1) + 1 : 8 * (j - 1) + 2, 15 + 12 * (tag_location - 1) + 1: 15 + 12 * (tag_location - 1) + 3) = ...
            Jp_1 * R_b_c * R_b_w';
        p2_world_coordinate = feature_history(12 * (tag_location - 1) + 4 : 12 * (tag_location - 1) + 6, 1);
        camera_coordinate(:, 4 * j - 2) = R_b_c * (R_b_w' * p2_world_coordinate - R_b_w' * T_b_w) - T_b_c;
        X2 = camera_coordinate(1, 4 * j - 2);
        Y2 = camera_coordinate(2, 4 * j - 2);
        Z2 = camera_coordinate(3, 4 * j - 2);
        predicted_measurement(:, 4 * j - 2) = [X2/Z2; Y2/Z2];
        Jp_2 = [1, 0, -X2/Z2; 0, 1, -Y2/Z2] / Z2;
        C(8 * (j - 1) + 3 : 8 * (j - 1) + 4, 1 : 3) = -Jp_2 * R_b_c * R_b_w';
        C(8 * (j - 1) + 3 : 8 * (j - 1) + 4, 7 : 9) = Jp_2 * R_b_c * R_b_w' * skew(p2_world_coordinate - T_b_w);
        C(8 * (j - 1) + 3 : 8 * (j - 1) + 4, 15 + 12 * (tag_location - 1) + 4: 15 + 12 * (tag_location - 1) + 6) = ...
            Jp_2 * R_b_c * R_b_w';
        p3_world_coordinate = feature_history(12 * (tag_location - 1) + 7 : 12 * (tag_location - 1) + 9, 1);
        camera_coordinate(:, 4 * j - 1) = R_b_c * (R_b_w' * p3_world_coordinate - R_b_w' * T_b_w) - T_b_c;
        X3 = camera_coordinate(1, 4 * j - 1);
        Y3 = camera_coordinate(2, 4 * j - 1);
        Z3 = camera_coordinate(3, 4 * j - 1);
        predicted_measurement(:, 4 * j - 1) = [X3/Z3; Y3/Z3];
        Jp_3 = [1, 0, -X3/Z3; 0, 1, -Y3/Z3] / Z3;
        C(8 * (j - 1) + 5 : 8 * (j - 1) + 6, 1 : 3) = -Jp_3 * R_b_c * R_b_w';
        C(8 * (j - 1) + 5 : 8 * (j - 1) + 6, 7 : 9) = Jp_3 * R_b_c * R_b_w' * skew(p3_world_coordinate - T_b_w);
        C(8 * (j - 1) + 5 : 8 * (j - 1) + 6, 15 + 12 * (tag_location - 1) + 7: 15 + 12 * (tag_location - 1) + 9) = ...
            Jp_3 * R_b_c * R_b_w';
        p4_world_coordinate = feature_history(12 * (tag_location - 1) + 10 : 12 * (tag_location - 1) + 12, 1);
        camera_coordinate(:, 4 * j) = R_b_c * (R_b_w' * p4_world_coordinate - R_b_w' * T_b_w) - T_b_c;
        X4 = camera_coordinate(1, 4 * j);
        Y4 = camera_coordinate(2, 4 * j);
        Z4 = camera_coordinate(3, 4 * j);
        predicted_measurement(:, 4 * j) = [X4/Z4; Y4/Z4];
        Jp_4 = [1, 0, -X4/Z4; 0, 1, -Y4/Z4] / Z4;
        C(8 * (j - 1) + 7 : 8 * (j - 1) + 8, 1 : 3) = -Jp_4 * R_b_c * R_b_w';
        C(8 * (j - 1) + 7 : 8 * (j - 1) + 8, 7 : 9) = Jp_4 * R_b_c * R_b_w' * skew(p4_world_coordinate - T_b_w);
        C(8 * (j - 1) + 7 : 8 * (j - 1) + 8, 15 + 12 * (tag_location - 1) + 10: 15 + 12 * (tag_location - 1) + 12) = ...
            Jp_4 * R_b_c * R_b_w';
    end
    actual_measurement = K \ actual_measurement;
    Y = reshape(actual_measurement(1:2, :), 8 * tag_number_seen, 1);
    X_bar = reshape(predicted_measurement, 8 * tag_number_seen, 1);
    V = 0.01 * eye(8 * tag_number_seen);
    K = sigma_together * C' / (C * sigma_together * C' + V);
    mu_correction = K * (Y - X_bar);
    sigma_correction = (eye(size(mu_correction, 1)) - K * C) * sigma_together * (eye(size(mu_correction, 1)) - K * C)' + K * V * K';
    dp_correction = mu_correction(1:3);
    dv_correction = mu_correction(4:6);
    d_theta_correction = mu_correction(7:9);
    d_bias_acceleration_correction = mu_correction(10:12);
    d_bias_gyro_correction = mu_correction(13:15);
    k_correction = d_theta_correction / norm(d_theta_correction);
    theta_correction = norm(d_theta_correction);
    dq_correction = [cos(theta_correction / 2); k_correction(1) * sin(theta_correction / 2); k_correction(2) * sin(theta_correction / 2); k_correction(3) * sin(theta_correction / 2)];
    p_correction = p_prediction + dp_correction;
    v_correction = v_prediction + dv_correction;
    q_correction = quatmultiply(q_prediction, dq_correction);
%     q_correction = q_correction';
    bias_acceleration_correction = bias_acceleration_prediction + d_bias_acceleration_correction;
    bias_gyro_correction = bias_gyro_prediction + d_bias_gyro_correction;
    p_history = p_correction;
    v_history = v_correction;
    q_history = q_correction;
    bias_acceleration_history = bias_acceleration_correction;
    bias_gyro_history = bias_gyro_correction;
    feature_history = feature_history + mu_correction(16:end);
    mu = zeros(15,1);
    hat_theta = skew(0.5 * d_theta);
    G = [eye(6),     zeros(6,3),         zeros(6,6);
         zeros(3,6), eye(3) - hat_theta, zeros(3,6);
         zeros(6,6), zeros(6,3),         eye(6)];
    sigma = G * sigma_correction(1:15,1:15) * G';
    feature_sigma = sigma_correction(16:end, 16:end);
end
%% Output
X = [p_history; v_history; q_history];
Z = [p_history;q_history];
end
