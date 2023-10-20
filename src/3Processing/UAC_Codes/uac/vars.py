#!/usr/bin/python3
#

EPSILON = 100
EPSILON_ABC = 0.025 * 0.01


# 21: LSPV. 23: LIPV. 25: RSPV. 27: RIPV

# inferior PV parameters
R1_23_27 = 0.06  # x_rad of inf PV
XSHIFT_23 = 0.07 + R1_23_27  # x_cen of LIPV
XSHIFT_27 = 1 - XSHIFT_23  # x_cen of RIPV
R2_23_27 = 2 * R1_23_27  # y_rad of inf PV
YSHIFT_23_27 = 0.63  # y_cen of inferior PV

R1_23_27_INNER = 0.02  # x_rad of inner inf PV
R2_23_27_INNER = 2 * R1_23_27_INNER  # y_rad of inner inf PV

LS_INF_PV_1 = XSHIFT_23 - R1_23_27
LS_INF_PV_2 = XSHIFT_23 + R1_23_27
LS_INF_PV_3 = XSHIFT_27 - R1_23_27
LS_INF_PV_4 = XSHIFT_27 + R1_23_27

PA_INF_PV_1 = YSHIFT_23_27 - R2_23_27
PA_INF_PV_2 = YSHIFT_23_27 + R2_23_27


# superior PV parameters
R1_21_25 = 0.115  # x_rad of sup PV
R2_21_25 = 0.22  # y_rad of sup PV

XSHIFT_21 = R1_21_25  # x_cen of LSPV
XSHIFT_25 = 1 - R1_21_25  # x_cen of RSPV
R1_21_25_INNER = 0.02  # x_rad of inner sup PV
R2_21_25_INNER = 2 * R1_21_25_INNER
YSHIFT_21_25 = 1 - R2_21_25 / 2  # y_cen of inner circles for superior PV

ANT_RSPV_LS = XSHIFT_25 - R1_21_25

# LAA parameters
R1_LAA = (R1_21_25 * 2) * 0.60 / 2  # x_rad of LAA
R2_LAA = 0.15  # y_rad of LAA
XSHIFT_LAA = R1_21_25 * 2 - R1_LAA  # x_cen of LAA
YSHIFT_LAA = 0.12 + R2_LAA  # y_cen of LAA

PA_ANT_LAA = YSHIFT_LAA + R2_LAA
PA_ANT_LAA_BOT = YSHIFT_LAA - R2_LAA

LS_ANT_LAA_1 = XSHIFT_LAA - R1_LAA
LS_ANT_LAA_2 = XSHIFT_LAA + R1_LAA
