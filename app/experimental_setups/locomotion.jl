using Voxcraft

biobot_size = (8, 8, 7)
min_cell_percentage = 0.3

EnableTemp()
BaseTemp(25)
TempAmp(39.4714242553)  # 50% volumetric change with temp_base=25: (1+0.01*(39.4714242553-25))**3-1=0.5
frequency = 5
TempPeriod(1/frequency)

init_time = 0.4
SimTime(3+init_time)

GravityAcc(-9.81)
LatticeDim(0.01)

EnableExpansion()
