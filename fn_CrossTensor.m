%Inputs: q: quaternion vector
%        flag:0-cross;1-star
function [CrossTensor] = fn_CrossTensor(q,flag)
    CrossTensor = zeros(4,4);
    q_v = q(1:3);
    q_0 = q(4);    
    cross_q_v = fn_VectorToSkewSymmetricTensor(q_v);
    switch (flag)
        case 0
            CrossTensor = [-cross_q_v + q_0*eye(3,3), q_v;-q_v',q_0];
        case 1
            CrossTensor = [cross_q_v + q_0*eye(3,3), q_v;-q_v',q_0];
        otherwise
            
    end
end
            