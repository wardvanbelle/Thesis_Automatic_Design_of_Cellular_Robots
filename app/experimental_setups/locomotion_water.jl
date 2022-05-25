using Voxcraft

biobot_size = (8, 8, 7)
min_cell_percentage = 0.3

EnableTemp()
BaseTemp(25)
TempAmp(39.4714242553)  # 50% volumetric change with temp_base=25: (1+0.01*(39.4714242553-25))**3-1=0.5
frequency = 4
TempPeriod(1/frequency)

init_time = 1
SimTime(10+init_time)

GravityAcc(-0.1)
LatticeDim(0.05)

EnableExpansion()

FluidEnv = true
RhoFluid = 1000 # water density
CDrag = 1.5 # fluid drag associated to a triangular facet
AggregateDragCoef = 0.5 * CDrag * RhoFluid # aggregate drag coefficient