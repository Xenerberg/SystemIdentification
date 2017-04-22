# SystemIdentification
Identification of Inertia/Geometry parameters of tumbling satellite

In the robotic missions of On-Orbit servicing (OOS) and Active debris removal (ADR), there is uncertainty in the dynamic/kinematic motion  parameters. Obviously, for the latter case, the uncertainty could be high enough to be unknown. In either case, it becomes important to robustly estimate/observe these parameters before any robotic commands are issued to initiate capturing. 

In this project, an Extended Kalman Filter was designed to perform this estimation. Specifically, two important motion parameters are observed. The first is the inertia of the tumbling body in unconstrained coordinates. Secondly, the pose of the observed frame on the body with respect to its center of mass is equally important for the dynamic model. The main purpose of this work is to establish a set of metrics that will allow comparison with alternate methods to achieve this overall objective.
