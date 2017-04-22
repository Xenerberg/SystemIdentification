function [psi] = fn_CreatePsi(q,p)
    psi = [-2*q(4)*fn_VectorToSkewSymmetricTensor(p) + 2*(fn_VectorToSkewSymmetricTensor(p)*fn_VectorToSkewSymmetricTensor(q(1:3)) - 2*fn_VectorToSkewSymmetricTensor(q(1:3))*fn_VectorToSkewSymmetricTensor(p)), 2*fn_VectorToSkewSymmetricTensor(q(1:3))*p];

end
