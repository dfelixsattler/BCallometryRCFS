# data-raw/biomass_coefs.R
#
# Builds the biomass_coefs package dataset.
# Sources:
#   Lambert M-C, Ung C-H, Raulier F (2005). Canadian national tree aboveground
#   biomass equations. Canadian Journal of Forest Research, 35(8), 1996-2018.
#
#   Ung C-H, Bernier P, Guo X-J (2008). Canadian national biomass equations:
#   new parameter estimates that include British Columbia data. Canadian Journal
#   of Forest Research, 38(5), 1123-1132.
#
# Coefficients transcribed from Table 3 (Lambert 2005) and the Ung 2008
# update. Coefficients were cross-verified against the original FAIBBase
# biomassCalculator implementation (github.com/bcgov/FAIBBase).
#
# Run this script once to regenerate data/biomass_coefs.rda.
# Do NOT edit data/biomass_coefs.rda directly.

# ---- helpers ----------------------------------------------------------------

# DBH-only equation:   B_component = a1 * DBH^a2          (a3 = NA)
# DBH+height equation: B_component = a1 * DBH^a2 * H^a3   (a3 = 0 means no
#                                                            height effect)

.d <- function(sp, ps,
               w1, w2, bk1, bk2, br1, br2, f1, f2) {
  data.frame(
    species         = sp,
    paper_source    = ps,
    height_included = FALSE,
    component       = c("wood", "bark", "branches", "foliage"),
    a1              = c(w1,  bk1, br1, f1),
    a2              = c(w2,  bk2, br2, f2),
    a3              = NA_real_,
    stringsAsFactors = FALSE
  )
}

.h <- function(sp, ps,
               w1, w2, w3, bk1, bk2, bk3, br1, br2, br3, f1, f2, f3) {
  data.frame(
    species         = sp,
    paper_source    = ps,
    height_included = TRUE,
    component       = c("wood", "bark", "branches", "foliage"),
    a1              = c(w1,  bk1, br1, f1),
    a2              = c(w2,  bk2, br2, f2),
    a3              = c(w3,  bk3, br3, f3),
    stringsAsFactors = FALSE
  )
}

L <- "Lambert2005"
U <- "Ung2008"

# ---- DBH-only entries (height_included = FALSE) -----------------------------
# Most species: Lambert2005 only (no Ung2008 branch in original source).
# Species with both L and U: black spruce, lodgepole pine, trembling aspen,
# white birch, white spruce, hardwood, softwood.

dbh_coefs <- rbind(
  .d("alpine fir",          L, 0.0528,2.4309, 0.0108,2.3876, 0.0121,2.3519, 0.0251,2.0389),
  .d("balsam fir",          L, 0.0534,2.4030, 0.0115,2.3484, 0.0070,2.5406, 0.0840,1.6695),
  .d("balsam poplar",       L, 0.0510,2.4529, 0.0297,2.1131, 0.0120,2.4165, 0.0276,1.6215),
  .d("basswood",            L, 0.0562,2.4102, 0.0302,2.0976, 0.0230,2.2382, 0.0288,1.6378),
  .d("black ash",           L, 0.0941,2.3491, 0.0323,2.0761, 0.0448,1.9771, 0.0538,1.3584),
  .d("black cherry",        L, 0.3743,1.9406, 0.0679,1.8377, 0.0796,2.0103, 0.0840,1.2319),
  .d("black spruce",        L, 0.0477,2.5147, 0.0153,2.2429, 0.0278,2.0839, 0.1648,1.4143),
  .d("black spruce",        U, 0.0494,2.5025, 0.0148,2.2494, 0.0291,2.0751, 0.1631,1.4222),
  .d("douglas-fir",         L, 0.0204,2.6974, 0.0069,2.5462, 0.0404,2.1388, 0.1233,1.6636),
  .d("engelmann spruce",    L, 0.0223,2.7169, 0.0118,2.2733, 0.0336,2.2123, 0.0683,1.8022),
  .d("eastern hemlock",     L, 0.0619,2.3821, 0.0139,2.3282, 0.0217,2.2653, 0.0776,1.6995),
  .d("eastern redcedar",    L, 0.1277,1.9778, 0.0377,1.6064, 0.0254,2.2884, 0.0550,1.8656),
  .d("eastern white-cedar", L, 0.0654,2.2121, 0.0114,2.1432, 0.0335,1.9367, 0.0499,1.7278),
  .d("eastern white pine",  L, 0.0997,2.2709, 0.0192,2.2038, 0.0056,2.6011, 0.0284,1.9375),
  .d("grey birch",          L, 0.0720,2.3885, 0.0168,2.2569, 0.0088,2.5689, 0.0099,1.8985),
  .d("hickory",             L, 0.2116,2.2013, 0.0365,2.1133, 0.0087,2.8927, 0.0173,1.9830),
  .d("hop-hornbeam",        L, 0.1929,1.9672, 0.0671,1.5911, 0.0278,2.1336, 0.0293,1.9502),
  .d("jack pine",           L, 0.0804,2.4041, 0.0184,2.0703, 0.0079,2.4155, 0.0389,1.7290),
  .d("largetooth aspen",    L, 0.0959,2.3430, 0.0308,2.2240, 0.0047,2.6530, 0.0080,2.0149),
  .d("lodgepole pine",      L, 0.0475,2.5437, 0.0186,2.0807, 0.0198,2.1287, 0.0432,1.7166),
  .d("lodgepole pine",      U, 0.0323,2.6825, 0.0144,2.1768, 0.0209,2.1772, 0.0584,1.6432),
  .d("pacific silver fir",  L, 0.0424,2.4289, 0.0057,2.4786, 0.0322,2.1313, 0.0645,1.9400),
  .d("red alder",           L, 0.0460,2.4312, 0.0074,2.4442, 0.0086,2.7326, 0.0114,2.0860),
  .d("black cottonwood",    L, 0.0460,2.4312, 0.0074,2.4442, 0.0086,2.7326, 0.0114,2.0860),
  .d("red ash",             L, 0.1571,2.1817, 0.0416,2.0509, 0.0177,2.3370, 0.1041,1.2185),
  .d("red maple",           L, 0.1014,2.3448, 0.0291,2.0893, 0.0175,2.4846, 0.0515,1.5198),
  .d("red oak",             L, 0.1754,2.1616, 0.0381,2.0991, 0.0085,2.7790, 0.0373,1.6740),
  .d("red pine",            L, 0.0564,2.4465, 0.0188,2.0527, 0.0033,2.7515, 0.0212,2.0690),
  .d("red spruce",          L, 0.0989,2.2814, 0.0220,2.0908, 0.0005,3.2750, 0.0066,2.4213),
  .d("sitka spruce",        L, 0.0302,2.5776, 0.0066,2.4433, 0.0739,1.8342, 0.0157,2.3113),
  .d("subalpine fir",       L, 0.0250,2.6378, 0.0061,2.5375, 0.0178,2.4255, 0.0416,2.0130),
  .d("silver maple",        L, 0.2324,2.1000, 0.0278,2.0433, 0.0028,3.1020, 0.1430,1.2580),
  .d("sugar maple",         L, 0.1315,2.3129, 0.0631,1.9241, 0.0330,2.3741, 0.0393,1.6930),
  .d("tamarack larch",      L, 0.0625,2.4475, 0.0174,2.1109, 0.0196,2.2652, 0.0801,1.4875),
  .d("trembling aspen",     L, 0.0605,2.4750, 0.0168,2.3949, 0.0080,2.5214, 0.0261,1.6304),
  .d("trembling aspen",     U, 0.0608,2.4735, 0.0159,2.4123, 0.0082,2.5139, 0.0235,1.6656),
  .d("western hemlock",     L, 0.0141,2.8668, 0.0025,2.8062, 0.0703,1.9547, 0.1676,1.4339),
  .d("western redcedar",    L, 0.0111,2.8027, 0.0003,3.2721, 0.1158,1.7196, 0.1233,1.5152),
  .d("white ash",           L, 0.1861,2.1665, 0.0406,1.9946, 0.0461,2.2291, 0.1106,1.2277),
  .d("white birch",         L, 0.0593,2.5026, 0.0135,2.4053, 0.0135,2.5532, 0.0546,1.6351),
  .d("white birch",         U, 0.0604,2.4959, 0.0140,2.3923, 0.0147,2.5227, 0.0591,1.6036),
  .d("white elm",           L, 0.0402,2.5804, 0.0073,2.4859, 0.0401,2.1826, 0.0750,1.3436),
  .d("white oak",           L, 0.0762,2.3335, 0.0338,1.9845, 0.0113,2.6211, 0.0188,1.7881),
  .d("white spruce",        L, 0.0359,2.5775, 0.0116,2.3022, 0.0283,2.0823, 0.1601,1.4670),
  .d("white spruce",        U, 0.0334,2.5980, 0.0114,2.3057, 0.0302,2.0927, 0.1515,1.5012),
  .d("yellow birch",        L, 0.1932,2.1569, 0.0192,2.2475, 0.0305,2.4044, 0.1119,1.3973),
  .d("hardwood",            L, 0.0871,2.3702, 0.0241,2.1969, 0.0167,2.4807, 0.0390,1.6229),
  .d("hardwood",            U, 0.0864,2.3715, 0.0226,2.2151, 0.0186,2.4462, 0.0385,1.6255),
  .d("softwood",            L, 0.0648,2.3927, 0.0162,2.1959, 0.0156,2.2916, 0.0861,1.6261),
  .d("softwood",            U, 0.0564,2.4347, 0.0153,2.2110, 0.0194,2.2408, 0.0935,1.6106)
)

# ---- DBH+height entries (height_included = TRUE) ----------------------------
# a3 = 0 means the component has no height effect (height^0 = 1).

ht_coefs <- rbind(
  # alpine fir — Lambert2005 only
  .h("alpine fir",L, 0.0268,1.7579,0.9871, 0.0009,1.4460,1.8839, 0.0470,2.9288,-1.1588, 0.0551,1.7585,0),
  # balsam fir — Lambert2005 only
  .h("balsam fir",L, 0.0294,1.8357,0.8640, 0.0053,2.0876,0.5842, 0.0117,3.5097,-1.3006, 0.1245,2.5230,-1.1230),
  # balsam poplar — Lambert2005 only
  .h("balsam poplar",L, 0.0117,1.7757,1.2555, 0.0180,1.8131,0.5144, 0.0112,3.0861,-0.7164, 0.0617,1.8615,-0.5375),
  # basswood — Lambert2005 only
  .h("basswood",L, 0.0168,1.9844,0.8989, 0.0057,1.5881,1.1472, 0.0039,2.0084,0.8588, 0.0147,1.8300,0),
  # beech — Lambert2005 only (no DBH-only equations available)
  .h("beech",L, 0.0432,2.0378,0.7000, 0.0049,1.9057,0.6770, 0.0355,2.3749,0, 0.0452,1.5567,0),
  # black ash — Lambert2005 only
  .h("black ash",L, 0.0306,2.1836,0.5740, 0.0897,2.2634,-0.5670, 0.0994,2.1630,-0.4809, 0.0124,1.0325,0.8747),
  # black cherry — Lambert2005 only
  .h("black cherry",L, 0.0181,1.7013,1.3057, 0.0101,1.5956,0.9190, 0.0005,2.8004,0.8603, 0.1976,1.4421,-0.5264),
  # black spruce — both
  .h("black spruce",L, 0.0309,1.7527,1.0014, 0.0115,1.7405,0.6589, 0.0380,3.2558,-1.4218, 0.2048,2.5754,-1.3704),
  .h("black spruce",U, 0.0335,1.7389,0.9835, 0.0132,1.7657,0.5775, 0.0405,3.1917,-1.3674, 0.2078,2.5517,-1.3453),
  # douglas-fir — Lambert2005 only
  .h("douglas-fir",L, 0.0191,1.5365,1.3634, 0.0083,2.4811,0, 0.0351,2.2421,0, 0.0718,2.2935,-0.4744),
  # engelmann spruce — Lambert2005 only
  .h("engelmann spruce",L, 0.0133,1.3303,1.6877, 0.0086,1.6216,0.8192, 0.0428,2.7965,-0.7328, 0.0854,2.4388,-0.7630),
  # eastern hemlock — Lambert2005 only
  .h("eastern hemlock",L, 0.0257,1.9277,0.8576, 0.0118,1.9893,0.4700, 0.0215,2.6553,-0.4682, 0.1471,2.0108,-0.6080),
  # eastern redcedar — Lambert2005 only
  .h("eastern redcedar",L, 0.0520,1.7731,0.7054, 0.0283,1.7079,0, 0.0219,2.3585,0, 0.2575,2.5136,-1.5565),
  # eastern white-cedar — Lambert2005 only
  .h("eastern white-cedar",L, 0.0295,1.7026,0.9428, 0.0076,1.7861,0.6132, 0.0501,2.5165,-0.8774, 0.0813,2.2180,-0.7907),
  # eastern white pine — Lambert2005 only
  .h("eastern white pine",L, 0.0170,1.7779,1.1370, 0.0069,1.6589,0.9582, 0.0184,3.1968,-1.0876, 0.0584,2.2389,-0.5968),
  # grey birch — Lambert2005 only
  .h("grey birch",L, 0.0295,1.9064,0.9139, 0.0148,1.8433,0.5021, 0.0150,3.0347,-0.7629, 0.0455,2.6447,-1.4955),
  # hickory — Lambert2005 only
  .h("hickory",L, 0.0139,1.5913,1.5080, 0.0081,1.4943,1.1324, 0.0050,3.0463,0, 0.0121,2.0865,0),
  # hop-hornbeam — Lambert2005 only
  .h("hop-hornbeam",L, 0.0083,1.6534,1.7479, 0.0012,1.1486,2.2903, 0.0009,1.9152,1.7769, 0.0247,2.0056,0),
  # jack pine — Lambert2005 only
  .h("jack pine",L, 0.0199,1.6883,1.2456, 0.0141,1.5994,0.5957, 0.0185,3.0584,-0.9816, 0.0325,1.7879,0),
  # largetooth aspen — Lambert2005 only
  .h("largetooth aspen",L, 0.0128,2.0633,0.9516, 0.0240,2.3055,0, 0.0131,3.1274,-0.8379, 0.0382,2.1673,-0.6842),
  # lodgepole pine — both
  .h("lodgepole pine",L, 0.0202,1.7179,1.2078, 0.0099,1.6049,0.7456, 0.0440,3.7190,-2.0399, 0.0785,2.5377,-1.1213),
  .h("lodgepole pine",U, 0.0239,1.6827,1.1878, 0.0117,1.6398,0.6524, 0.0285,3.3764,-1.4395, 0.0769,2.6834,-1.2484),
  # pacific silver fir — Lambert2005 only
  .h("pacific silver fir",L, 0.0315,1.8297,0.8056, 0.0067,2.6970,-0.3105, 0.0420,2.0313,0, 0.0453,2.4867,-0.4982),
  # red alder / black cottonwood — Lambert2005 only (same coefficients)
  .h("red alder",L, 0.0051,1.0697,2.2748, 0.0009,1.3061,2.0109, 0.0131,2.5760,0, 0.0224,1.8368,0),
  .h("black cottonwood",L, 0.0051,1.0697,2.2748, 0.0009,1.3061,2.0109, 0.0131,2.5760,0, 0.0224,1.8368,0),
  # red ash — Lambert2005 only
  .h("red ash",L, 0.0224,1.7845,1.0600, 0.0219,1.4190,0.8963, 0.0176,2.3313,0, 0.0761,1.3077,0),
  # red maple — Lambert2005 only
  .h("red maple",L, 0.0315,2.0342,0.7485, 0.0283,2.0907,0, 0.0225,2.4106,0, 0.0571,1.4898,0),
  # red oak — Lambert2005 only
  .h("red oak",L, 0.0285,1.8501,1.0204, 0.0326,1.8100,0.4153, 0.0013,3.0637,0.3153, 0.0582,1.5438,0),
  # red pine — Lambert2005 only
  .h("red pine",L, 0.0106,1.7725,1.3285, 0.0277,1.5192,0.4645, 0.0125,3.3865,-1.1939, 0.0731,2.3439,-0.7378),
  # red spruce — Lambert2005 only
  .h("red spruce",L, 0.0143,1.6441,1.4065, 0.0274,2.0188,0, 0.0005,3.3136,0, 0.0106,2.2709,0),
  # sitka spruce — Lambert2005 only
  .h("sitka spruce",L, 0.0237,2.5813,0.0822, 0.0045,1.2275,1.5190, 0.0498,1.9671,0, 0.0140,3.1305,-0.9070),
  # subalpine fir — Lambert2005 only
  .h("subalpine fir",L, 0.0220,1.6469,1.1714, 0.0061,1.8603,0.7693, 0.0265,3.6747,-1.5958, 0.0509,2.9909,-1.2271),
  # silver maple — Lambert2005 only
  .h("silver maple",L, 0.0274,1.7126,1.1086, 0.0123,1.8250,0.5010, 0.0543,3.7343,-1.6497, 6.6808,2.1092,-2.1697),
  # sugar maple — Lambert2005 only
  .h("sugar maple",L, 0.0301,2.0313,0.8171, 0.0103,1.7111,0.8509, 0.0661,2.5940,-0.4933, 2.5019,2.4527,-2.3008),
  # tamarack larch — Lambert2005 only
  .h("tamarack larch",L, 0.0276,1.6724,1.1443, 0.0120,1.7059,0.5811, 0.0336,3.1335,-1.1559, 0.1324,2.1140,-0.8781),
  # trembling aspen — both
  .h("trembling aspen",L, 0.0142,1.9389,1.0572, 0.0063,2.0819,0.6617, 0.0137,2.9270,-0.6221, 0.0270,1.6183,0),
  .h("trembling aspen",U, 0.0143,1.9369,1.0579, 0.0063,2.0744,0.6691, 0.0150,2.9068,-0.6306, 0.0284,1.6020,0),
  # western hemlock — Lambert2005 only
  .h("western hemlock",L, 0.0113,1.9332,1.1125, 0.0019,2.3356,0.6371, 0.0609,2.0021,0, 0.2656,2.0107,-0.7963),
  # western redcedar — Lambert2005 only
  .h("western redcedar",L, 0.0188,1.3376,1.5293, 0.0002,2.4369,1.1315, 0.0611,1.9208,0, 0.1097,1.5530,0),
  # white ash — Lambert2005 only
  .h("white ash",L, 0.0224,1.7438,1.1899, 0.0126,1.6456,0.7893, 0.0354,2.3046,0, 0.0195,1.0509,0.7836),
  # white birch — both
  .h("white birch",L, 0.0338,2.0702,0.6876, 0.0080,1.9754,0.6659, 0.0257,3.1754,-0.9417, 0.1415,2.3074,-1.1189),
  .h("white birch",U, 0.0333,2.0794,0.6811, 0.0079,1.9905,0.6553, 0.0253,3.1518,-0.9083, 0.1361,2.2978,-1.0934),
  # white elm — Lambert2005 only
  .h("white elm",L, 0.0207,2.2276,0.6488, 0.0078,2.4540,0, 0.0393,2.1880,0, 0.0516,1.4511,0),
  # white oak — Lambert2005 only
  .h("white oak",L, 0.0442,1.6818,1.0310, 0.0308,1.7479,0.3504, 0.0022,2.0165,1.3953, 0.0053,1.2822,1.1323),
  # white spruce — both
  .h("white spruce",L, 0.0265,1.7952,0.9733, 0.0124,1.6962,0.6489, 0.0325,2.8573,-0.9127, 0.2020,2.3802,-1.1103),
  .h("white spruce",U, 0.0252,1.7819,1.0022, 0.0096,1.6901,0.7393, 0.0322,2.8961,-0.9203, 0.1832,2.4144,-1.0948),
  # yellow birch — Lambert2005 only
  .h("yellow birch",L, 0.0259,1.9044,0.9715, 0.0069,2.0834,0.5371, 0.0325,2.3851,0, 0.1683,1.2764,0),
  # hardwood generic — both
  .h("hardwood",L, 0.0359,2.0263,0.6987, 0.0094,1.8677,0.6985, 0.0433,2.6817,-0.5731, 0.0859,1.8485,-0.5383),
  .h("hardwood",U, 0.0353,2.0249,0.7048, 0.0090,1.8677,0.7144, 0.0448,2.6855,-0.5911, 0.0869,1.8541,-0.5491),
  # softwood generic — both
  .h("softwood",L, 0.0284,1.6894,1.0857, 0.0100,1.8463,0.5616, 0.0301,3.0038,-1.0520, 0.1554,2.4021,-1.1043),
  .h("softwood",U, 0.0276,1.6868,1.0953, 0.0101,1.8486,0.5525, 0.0313,2.9974,-1.0383, 0.1379,2.3981,-1.0418),
  # generic fallback (used when species not recognised, height_included = TRUE)
  .h("generic",L, 0.0348,1.9235,0.7829, 0.0139,1.5429,0.8189, 0.0346,2.6706,-0.6033, 0.1822,2.2864,-1.1203),
  .h("generic",U, 0.0283,1.8298,0.9546, 0.0120,1.6378,0.7746, 0.0338,2.6624,-0.5743, 0.1699,2.3289,-1.1316)
)

# ---- combine and save -------------------------------------------------------
biomass_coefs <- rbind(dbh_coefs, ht_coefs)

cat(sprintf("biomass_coefs: %d rows, %d species, %d paper sources\n",
            nrow(biomass_coefs),
            length(unique(biomass_coefs$species)),
            length(unique(biomass_coefs$paper_source))))

usethis::use_data(biomass_coefs, overwrite = TRUE)
