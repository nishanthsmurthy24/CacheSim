#!/bin/bash -l
#
#
#
# Traces:
#    403.gcc-16B_4T
#    410.bwaves-1963B_4T
#    429.mcf-184B_4T
#    429.mcf-192B_4T
#    429.mcf-51B_4T
#    433.milc-127B_4T
#    436.cactusADM-1804B_4T
#    437.leslie3d-149B_4T
#    437.leslie3d-232B_4T
#    437.leslie3d-265B_4T
#    437.leslie3d-271B_4T
#    445.gobmk-30B_4T
#    445.gobmk-36B_4T
#    450.soplex-247B_4T
#    459.GemsFDTD-1211B_4T
#    459.GemsFDTD-1418B_4T
#    459.GemsFDTD-765B_4T
#    462.libquantum-1343B_4T
#    462.libquantum-714B_4T
#    470.lbm-1274B_4T
#    471.omnetpp-188B_4T
#    473.astar-359B_4T
#    481.wrf-1254B_4T
#    481.wrf-196B_4T
#    481.wrf-816B_4T
#    482.sphinx3-1297B_4T
#    482.sphinx3-1395B_4T
#    483.xalancbmk-127B_4T
#    602.gcc_s-1850B_4T
#    602.gcc_s-734B_4T
#    603.bwaves_s-2931B_4T
#    605.mcf_s-994B_4T
#    607.cactuBSSN_s-2421B_4T
#    619.lbm_s-2677B_4T
#    619.lbm_s-3766B_4T
#    619.lbm_s-4268B_4T
#    620.omnetpp_s-141B_4T
#    620.omnetpp_s-874B_4T
#    621.wrf_s-6673B_4T
#    621.wrf_s-8065B_4T
#    623.xalancbmk_s-10B_4T
#    623.xalancbmk_s-592B_4T
#    627.cam4_s-573B_4T
#    628.pop2_s-17B_4T
#    649.fotonik3d_s-10881B_4T
#    654.roms_s-1007B_4T
#    ligra_BC.com-lj.ungraph.gcc_6.3.0_O3.drop_26750M.length_250M_4T
#    ligra_BC.com-lj.ungraph.gcc_6.3.0_O3.drop_3500M.length_250M_4T
#    ligra_BC.com-lj.ungraph.gcc_6.3.0_O3.drop_500M.length_250M_4T
#    ligra_BellmanFord.com-lj.ungraph.gcc_6.3.0_O3.drop_1750M.length_250M_4T
#    ligra_BellmanFord.com-lj.ungraph.gcc_6.3.0_O3.drop_4000M.length_250M_4T
#    ligra_BellmanFord.com-lj.ungraph.gcc_6.3.0_O3.drop_7500M.length_250M_4T
#    ligra_BFS-Bitvector.com-lj.ungraph.gcc_6.3.0_O3.drop_23000M.length_250M_4T
#    ligra_BFS-Bitvector.com-lj.ungraph.gcc_6.3.0_O3.drop_2500M.length_250M_4T
#    ligra_BFSCC.com-lj.ungraph.gcc_6.3.0_O3.drop_22000M.length_250M_4T
#    ligra_BFSCC.com-lj.ungraph.gcc_6.3.0_O3.drop_3500M.length_250M_4T
#    ligra_BFSCC.com-lj.ungraph.gcc_6.3.0_O3.drop_5000M.length_250M_4T
#    ligra_BFSCC.com-lj.ungraph.gcc_6.3.0_O3.drop_750M.length_250M_4T
#    ligra_BFS.com-lj.ungraph.gcc_6.3.0_O3.drop_21500M.length_250M_4T
#    ligra_BFS.com-lj.ungraph.gcc_6.3.0_O3.drop_3500M.length_250M_4T
#    ligra_BFS.com-lj.ungraph.gcc_6.3.0_O3.drop_5000M.length_250M_4T
#    ligra_BFS.com-lj.ungraph.gcc_6.3.0_O3.drop_500M.length_250M_4T
#    ligra_CF.com-lj.ungraph.gcc_6.3.0_O3.drop_184750M.length_250M_4T
#    ligra_CF.com-lj.ungraph.gcc_6.3.0_O3.drop_2500M.length_250M_4T
#    ligra_Components.com-lj.ungraph.gcc_6.3.0_O3.drop_22750M.length_250M_4T
#    ligra_Components.com-lj.ungraph.gcc_6.3.0_O3.drop_3500M.length_250M_4T
#    ligra_Components.com-lj.ungraph.gcc_6.3.0_O3.drop_750M.length_250M_4T
#    ligra_Components-Shortcut.com-lj.ungraph.gcc_6.3.0_O3.drop_22000M.length_250M_4T
#    ligra_Components-Shortcut.com-lj.ungraph.gcc_6.3.0_O3.drop_750M.length_250M_4T
#    ligra_MIS.com-lj.ungraph.gcc_6.3.0_O3.drop_21250M.length_250M_4T
#    ligra_MIS.com-lj.ungraph.gcc_6.3.0_O3.drop_3500M.length_250M_4T
#    ligra_MIS.com-lj.ungraph.gcc_6.3.0_O3.drop_750M.length_250M_4T
#    ligra_PageRank.com-lj.ungraph.gcc_6.3.0_O3.drop_21750M.length_250M_4T
#    ligra_PageRank.com-lj.ungraph.gcc_6.3.0_O3.drop_500M.length_250M_4T
#    ligra_PageRank.com-lj.ungraph.gcc_6.3.0_O3.drop_51000M.length_250M_4T
#    ligra_PageRank.com-lj.ungraph.gcc_6.3.0_O3.drop_60750M.length_250M_4T
#    ligra_PageRank.com-lj.ungraph.gcc_6.3.0_O3.drop_79500M.length_250M_4T
#    ligra_PageRankDelta.com-lj.ungraph.gcc_6.3.0_O3.drop_1250M.length_250M_4T
#    ligra_PageRankDelta.com-lj.ungraph.gcc_6.3.0_O3.drop_24000M.length_250M_4T
#    ligra_PageRankDelta.com-lj.ungraph.gcc_6.3.0_O3.drop_24500M.length_250M_4T
#    ligra_PageRankDelta.com-lj.ungraph.gcc_6.3.0_O3.drop_3500M.length_250M_4T
#    ligra_Radii.com-lj.ungraph.gcc_6.3.0_O3.drop_32000M.length_250M_4T
#    ligra_Radii.com-lj.ungraph.gcc_6.3.0_O3.drop_3500M.length_250M_4T
#    ligra_Triangle.com-lj.ungraph.gcc_6.3.0_O3.drop_25000M.length_250M_4T
#    ligra_Triangle.com-lj.ungraph.gcc_6.3.0_O3.drop_3500M.length_250M_4T
#    ligra_Triangle.com-lj.ungraph.gcc_6.3.0_O3.drop_750M.length_250M_4T
#    parsec_2.1.canneal.simlarge.prebuilt.drop_1250M.length_250M_4T
#    parsec_2.1.canneal.simlarge.prebuilt.drop_3000M.length_250M_4T
#    parsec_2.1.canneal.simlarge.prebuilt.drop_4500M.length_250M_4T
#    parsec_2.1.canneal.simlarge.prebuilt.drop_4750M.length_250M_4T
#    parsec_2.1.canneal.simlarge.prebuilt.drop_500M.length_250M_4T
#    parsec_2.1.facesim.simlarge.prebuilt.drop_1500M.length_250M_4T
#    parsec_2.1.fluidanimate.simlarge.prebuilt.drop_9500M.length_250M_4T
#    parsec_2.1.streamcluster.simlarge.prebuilt.drop_0M.length_250M_4T
#    parsec_2.1.streamcluster.simlarge.prebuilt.drop_14750M.length_250M_4T
#    parsec_2.1.raytrace.simlarge.prebuilt.drop_23500M.length_250M_4T
#    parsec_2.1.raytrace.simlarge.prebuilt.drop_23750M.length_250M_4T
#    cassandra_phase0_core0_4T
#    cassandra_phase0_core1_4T
#    cassandra_phase0_core2_4T
#    cassandra_phase0_core3_4T
#    cassandra_phase1_core0_4T
#    cassandra_phase1_core1_4T
#    cassandra_phase1_core2_4T
#    cassandra_phase1_core3_4T
#    cassandra_phase2_core1_4T
#    cassandra_phase2_core2_4T
#    cassandra_phase2_core3_4T
#    cassandra_phase3_core1_4T
#    cassandra_phase3_core3_4T
#    cassandra_phase4_core0_4T
#    cassandra_phase4_core2_4T
#    cassandra_phase4_core3_4T
#    cassandra_phase5_core0_4T
#    cassandra_phase5_core1_4T
#    cassandra_phase5_core2_4T
#    cassandra_phase5_core3_4T
#    cloud9_phase5_core2_4T
#    nutch_phase0_core0_4T
#    nutch_phase0_core1_4T
#    nutch_phase0_core2_4T
#    nutch_phase0_core3_4T
#    nutch_phase1_core0_4T
#    nutch_phase1_core1_4T
#    nutch_phase1_core2_4T
#    nutch_phase1_core3_4T
#    nutch_phase3_core0_4T
#    nutch_phase3_core1_4T
#    nutch_phase3_core2_4T
#    nutch_phase3_core3_4T
#    nutch_phase4_core0_4T
#    nutch_phase4_core1_4T
#    nutch_phase4_core2_4T
#    nutch_phase4_core3_4T
#    streaming_phase0_core1_4T
#    streaming_phase1_core0_4T
#    streaming_phase1_core1_4T
#    streaming_phase1_core3_4T
#    streaming_phase2_core0_4T
#    streaming_phase2_core1_4T
#    streaming_phase2_core2_4T
#    streaming_phase2_core3_4T
#    streaming_phase3_core0_4T
#    streaming_phase3_core1_4T
#    streaming_phase3_core3_4T
#    streaming_phase4_core0_4T
#    streaming_phase4_core1_4T
#    streaming_phase4_core3_4T
#    streaming_phase5_core0_4T
#    streaming_phase5_core1_4T
#    MP_mix_2
#    MP_mix_5
#    MP_mix_7
#    MP_mix_10
#    MP_mix_11
#    MP_mix_12
#    MP_mix_15
#    MP_mix_16
#    MP_mix_18
#    MP_mix_19
#    MP_mix_21
#    MP_mix_23
#    MP_mix_24
#    MP_mix_25
#    MP_mix_29
#    MP_mix_32
#    MP_mix_34
#    MP_mix_36
#    MP_mix_39
#    MP_mix_40
#    MP_mix_41
#    MP_mix_42
#    MP_mix_45
#    MP_mix_46
#    MP_mix_49
#    MP_mix_52
#    MP_mix_54
#    MP_mix_55
#    MP_mix_57
#    MP_mix_59
#    MP_mix_60
#    MP_mix_61
#    MP_mix_62
#    MP_mix_63
#    MP_mix_64
#    MP_mix_65
#    MP_mix_66
#    MP_mix_67
#    MP_mix_69
#    MP_mix_70
#    MP_mix_74
#    MP_mix_75
#    MP_mix_76
#    MP_mix_77
#    MP_mix_79
#    MP_mix_82
#    MP_mix_84
#    MP_mix_87
#    MP_mix_89
#    MP_mix_91
#    MP_mix_92
#    MP_mix_93
#    MP_mix_95
#    MP_mix_96
#    MP_mix_97
#    MP_mix_100
#    MP_mix_101
#    MP_mix_102
#    MP_mix_103
#    MP_mix_104
#    MP_mix_106
#    MP_mix_107
#    MP_mix_110
#    MP_mix_111
#    MP_mix_116
#    MP_mix_117
#    MP_mix_121
#    MP_mix_123
#    MP_mix_124
#    MP_mix_126
#    MP_mix_127
#    MP_mix_128
#    MP_mix_134
#    MP_mix_136
#    MP_mix_137
#    MP_mix_138
#    MP_mix_139
#    MP_mix_142
#    MP_mix_143
#    MP_mix_145
#    MP_mix_147
#    MP_mix_150
#    MP_mix_151
#    MP_mix_152
#    MP_mix_154
#    MP_mix_158
#    MP_mix_160
#    MP_mix_161
#    MP_mix_162
#    MP_mix_164
#    MP_mix_165
#    MP_mix_166
#    MP_mix_167
#    MP_mix_168
#    MP_mix_169
#    MP_mix_170
#    MP_mix_171
#    MP_mix_172
#    MP_mix_175
#    MP_mix_178
#    MP_mix_179
#    MP_mix_181
#    MP_mix_182
#    MP_mix_183
#    MP_mix_188
#    MP_mix_191
#    MP_mix_195
#    MP_mix_197
#    MP_mix_198
#    MP_mix_199
#    MP_mix_200
#    MP_mix_201
#    MP_mix_202
#    MP_mix_205
#    MP_mix_207
#    MP_mix_208
#    MP_mix_211
#    MP_mix_212
#    MP_mix_215
#    MP_mix_220
#    MP_mix_221
#    MP_mix_222
#    MP_mix_223
#    MP_mix_224
#    MP_mix_226
#    MP_mix_227
#    MP_mix_228
#    MP_mix_229
#    MP_mix_230
#    MP_mix_231
#    MP_mix_233
#    MP_mix_235
#    MP_mix_238
#    MP_mix_240
#    MP_mix_241
#    MP_mix_244
#    MP_mix_246
#    MP_mix_247
#    MP_mix_253
#    MP_mix_255
#    MP_mix_258
#    MP_mix_260
#    MP_mix_262
#    MP_mix_263
#    MP_mix_266
#    MP_mix_273
#    MP_mix_274
#    MP_mix_275
#    MP_mix_277
#    MP_mix_278
#
#
# Experiments:
#    nopref: --warmup_instructions=50000000 --simulation_instructions=150000000 --config=$(PYTHIA_HOME)/config/nopref.ini
#    spp: --warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=spp_dev2 --config=$(PYTHIA_HOME)/config/spp_dev2.ini
#    stride: --warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=stride --config=$(PYTHIA_HOME)/config/stride.ini
#    bingo: --warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=bingo --config=$(PYTHIA_HOME)/config/bingo.ini
#    mlop: --warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=mlop --config=$(PYTHIA_HOME)/config/mlop.ini
#    dspatch: --warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=dspatch --config=$(PYTHIA_HOME)/config/dspatch.ini
#    pythia: --warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=scooby --config=$(PYTHIA_HOME)/config/pythia.ini
#
#
#
#
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_2_spp -o MP_mix_2_spp.out -e MP_mix_2_spp.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=spp_dev2 --config=/home/sartorij/somas026/CacheSim/config/spp_dev2.ini  -traces /home/sartorij/somas026/CacheSim/traces/437.leslie3d-265B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/621.wrf_s-8065B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/429.mcf-192B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/481.wrf-816B.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_2_stride -o MP_mix_2_stride.out -e MP_mix_2_stride.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=stride --config=/home/sartorij/somas026/CacheSim/config/stride.ini  -traces /home/sartorij/somas026/CacheSim/traces/437.leslie3d-265B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/621.wrf_s-8065B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/429.mcf-192B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/481.wrf-816B.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_2_bingo -o MP_mix_2_bingo.out -e MP_mix_2_bingo.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=bingo --config=/home/sartorij/somas026/CacheSim/config/bingo.ini  -traces /home/sartorij/somas026/CacheSim/traces/437.leslie3d-265B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/621.wrf_s-8065B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/429.mcf-192B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/481.wrf-816B.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_2_mlop -o MP_mix_2_mlop.out -e MP_mix_2_mlop.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=mlop --config=/home/sartorij/somas026/CacheSim/config/mlop.ini  -traces /home/sartorij/somas026/CacheSim/traces/437.leslie3d-265B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/621.wrf_s-8065B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/429.mcf-192B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/481.wrf-816B.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_2_pythia -o MP_mix_2_pythia.out -e MP_mix_2_pythia.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=scooby --config=/home/sartorij/somas026/CacheSim/config/pythia.ini  -traces /home/sartorij/somas026/CacheSim/traces/437.leslie3d-265B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/621.wrf_s-8065B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/429.mcf-192B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/481.wrf-816B.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_5_nopref -o MP_mix_5_nopref.out -e MP_mix_5_nopref.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --config=/home/sartorij/somas026/CacheSim/config/nopref.ini  -traces /home/sartorij/somas026/CacheSim/traces/654.roms_s-1613B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/619.lbm_s-2676B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/459.GemsFDTD-1211B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/429.mcf-184B.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_5_spp -o MP_mix_5_spp.out -e MP_mix_5_spp.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=spp_dev2 --config=/home/sartorij/somas026/CacheSim/config/spp_dev2.ini  -traces /home/sartorij/somas026/CacheSim/traces/654.roms_s-1613B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/619.lbm_s-2676B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/459.GemsFDTD-1211B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/429.mcf-184B.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_5_stride -o MP_mix_5_stride.out -e MP_mix_5_stride.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=stride --config=/home/sartorij/somas026/CacheSim/config/stride.ini  -traces /home/sartorij/somas026/CacheSim/traces/654.roms_s-1613B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/619.lbm_s-2676B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/459.GemsFDTD-1211B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/429.mcf-184B.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_5_bingo -o MP_mix_5_bingo.out -e MP_mix_5_bingo.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=bingo --config=/home/sartorij/somas026/CacheSim/config/bingo.ini  -traces /home/sartorij/somas026/CacheSim/traces/654.roms_s-1613B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/619.lbm_s-2676B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/459.GemsFDTD-1211B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/429.mcf-184B.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_5_mlop -o MP_mix_5_mlop.out -e MP_mix_5_mlop.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=mlop --config=/home/sartorij/somas026/CacheSim/config/mlop.ini  -traces /home/sartorij/somas026/CacheSim/traces/654.roms_s-1613B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/619.lbm_s-2676B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/459.GemsFDTD-1211B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/429.mcf-184B.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_5_dspatch -o MP_mix_5_dspatch.out -e MP_mix_5_dspatch.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=dspatch --config=/home/sartorij/somas026/CacheSim/config/dspatch.ini  -traces /home/sartorij/somas026/CacheSim/traces/654.roms_s-1613B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/619.lbm_s-2676B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/459.GemsFDTD-1211B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/429.mcf-184B.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_5_pythia -o MP_mix_5_pythia.out -e MP_mix_5_pythia.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=scooby --config=/home/sartorij/somas026/CacheSim/config/pythia.ini  -traces /home/sartorij/somas026/CacheSim/traces/654.roms_s-1613B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/619.lbm_s-2676B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/459.GemsFDTD-1211B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/429.mcf-184B.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_10_spp -o MP_mix_10_spp.out -e MP_mix_10_spp.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=spp_dev2 --config=/home/sartorij/somas026/CacheSim/config/spp_dev2.ini  -traces /home/sartorij/somas026/CacheSim/traces/437.leslie3d-273B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/429.mcf-192B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/429.mcf-184B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/470.lbm-1274B.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_10_stride -o MP_mix_10_stride.out -e MP_mix_10_stride.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=stride --config=/home/sartorij/somas026/CacheSim/config/stride.ini  -traces /home/sartorij/somas026/CacheSim/traces/437.leslie3d-273B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/429.mcf-192B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/429.mcf-184B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/470.lbm-1274B.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_10_bingo -o MP_mix_10_bingo.out -e MP_mix_10_bingo.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=bingo --config=/home/sartorij/somas026/CacheSim/config/bingo.ini  -traces /home/sartorij/somas026/CacheSim/traces/437.leslie3d-273B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/429.mcf-192B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/429.mcf-184B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/470.lbm-1274B.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_10_mlop -o MP_mix_10_mlop.out -e MP_mix_10_mlop.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=mlop --config=/home/sartorij/somas026/CacheSim/config/mlop.ini  -traces /home/sartorij/somas026/CacheSim/traces/437.leslie3d-273B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/429.mcf-192B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/429.mcf-184B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/470.lbm-1274B.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_10_dspatch -o MP_mix_10_dspatch.out -e MP_mix_10_dspatch.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=dspatch --config=/home/sartorij/somas026/CacheSim/config/dspatch.ini  -traces /home/sartorij/somas026/CacheSim/traces/437.leslie3d-273B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/429.mcf-192B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/429.mcf-184B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/470.lbm-1274B.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_10_pythia -o MP_mix_10_pythia.out -e MP_mix_10_pythia.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=scooby --config=/home/sartorij/somas026/CacheSim/config/pythia.ini  -traces /home/sartorij/somas026/CacheSim/traces/437.leslie3d-273B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/429.mcf-192B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/429.mcf-184B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/470.lbm-1274B.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_15_nopref -o MP_mix_15_nopref.out -e MP_mix_15_nopref.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --config=/home/sartorij/somas026/CacheSim/config/nopref.ini  -traces /home/sartorij/somas026/CacheSim/traces/627.cam4_s-573B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/654.roms_s-1007B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/619.lbm_s-2676B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/429.mcf-192B.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_15_spp -o MP_mix_15_spp.out -e MP_mix_15_spp.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=spp_dev2 --config=/home/sartorij/somas026/CacheSim/config/spp_dev2.ini  -traces /home/sartorij/somas026/CacheSim/traces/627.cam4_s-573B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/654.roms_s-1007B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/619.lbm_s-2676B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/429.mcf-192B.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_15_stride -o MP_mix_15_stride.out -e MP_mix_15_stride.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=stride --config=/home/sartorij/somas026/CacheSim/config/stride.ini  -traces /home/sartorij/somas026/CacheSim/traces/627.cam4_s-573B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/654.roms_s-1007B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/619.lbm_s-2676B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/429.mcf-192B.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_15_bingo -o MP_mix_15_bingo.out -e MP_mix_15_bingo.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=bingo --config=/home/sartorij/somas026/CacheSim/config/bingo.ini  -traces /home/sartorij/somas026/CacheSim/traces/627.cam4_s-573B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/654.roms_s-1007B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/619.lbm_s-2676B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/429.mcf-192B.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_15_mlop -o MP_mix_15_mlop.out -e MP_mix_15_mlop.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=mlop --config=/home/sartorij/somas026/CacheSim/config/mlop.ini  -traces /home/sartorij/somas026/CacheSim/traces/627.cam4_s-573B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/654.roms_s-1007B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/619.lbm_s-2676B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/429.mcf-192B.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_15_dspatch -o MP_mix_15_dspatch.out -e MP_mix_15_dspatch.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=dspatch --config=/home/sartorij/somas026/CacheSim/config/dspatch.ini  -traces /home/sartorij/somas026/CacheSim/traces/627.cam4_s-573B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/654.roms_s-1007B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/619.lbm_s-2676B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/429.mcf-192B.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_15_pythia -o MP_mix_15_pythia.out -e MP_mix_15_pythia.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=scooby --config=/home/sartorij/somas026/CacheSim/config/pythia.ini  -traces /home/sartorij/somas026/CacheSim/traces/627.cam4_s-573B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/654.roms_s-1007B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/619.lbm_s-2676B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/429.mcf-192B.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_19_spp -o MP_mix_19_spp.out -e MP_mix_19_spp.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=spp_dev2 --config=/home/sartorij/somas026/CacheSim/config/spp_dev2.ini  -traces /home/sartorij/somas026/CacheSim/traces/473.astar-359B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/605.mcf_s-1644B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/621.wrf_s-6673B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/619.lbm_s-4268B.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_19_bingo -o MP_mix_19_bingo.out -e MP_mix_19_bingo.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=bingo --config=/home/sartorij/somas026/CacheSim/config/bingo.ini  -traces /home/sartorij/somas026/CacheSim/traces/473.astar-359B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/605.mcf_s-1644B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/621.wrf_s-6673B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/619.lbm_s-4268B.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_19_mlop -o MP_mix_19_mlop.out -e MP_mix_19_mlop.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=mlop --config=/home/sartorij/somas026/CacheSim/config/mlop.ini  -traces /home/sartorij/somas026/CacheSim/traces/473.astar-359B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/605.mcf_s-1644B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/621.wrf_s-6673B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/619.lbm_s-4268B.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_19_dspatch -o MP_mix_19_dspatch.out -e MP_mix_19_dspatch.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=dspatch --config=/home/sartorij/somas026/CacheSim/config/dspatch.ini  -traces /home/sartorij/somas026/CacheSim/traces/473.astar-359B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/605.mcf_s-1644B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/621.wrf_s-6673B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/619.lbm_s-4268B.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_24_spp -o MP_mix_24_spp.out -e MP_mix_24_spp.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=spp_dev2 --config=/home/sartorij/somas026/CacheSim/config/spp_dev2.ini  -traces /home/sartorij/somas026/CacheSim/traces/654.roms_s-1007B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/603.bwaves_s-1740B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/473.astar-359B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/481.wrf-196B.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_24_mlop -o MP_mix_24_mlop.out -e MP_mix_24_mlop.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=mlop --config=/home/sartorij/somas026/CacheSim/config/mlop.ini  -traces /home/sartorij/somas026/CacheSim/traces/654.roms_s-1007B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/603.bwaves_s-1740B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/473.astar-359B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/481.wrf-196B.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_24_pythia -o MP_mix_24_pythia.out -e MP_mix_24_pythia.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=scooby --config=/home/sartorij/somas026/CacheSim/config/pythia.ini  -traces /home/sartorij/somas026/CacheSim/traces/654.roms_s-1007B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/603.bwaves_s-1740B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/473.astar-359B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/481.wrf-196B.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_41_spp -o MP_mix_41_spp.out -e MP_mix_41_spp.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=spp_dev2 --config=/home/sartorij/somas026/CacheSim/config/spp_dev2.ini  -traces /home/sartorij/somas026/CacheSim/traces/429.mcf-184B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/429.mcf-184B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/437.leslie3d-232B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/437.leslie3d-273B.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_41_bingo -o MP_mix_41_bingo.out -e MP_mix_41_bingo.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=bingo --config=/home/sartorij/somas026/CacheSim/config/bingo.ini  -traces /home/sartorij/somas026/CacheSim/traces/429.mcf-184B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/429.mcf-184B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/437.leslie3d-232B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/437.leslie3d-273B.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_41_pythia -o MP_mix_41_pythia.out -e MP_mix_41_pythia.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=scooby --config=/home/sartorij/somas026/CacheSim/config/pythia.ini  -traces /home/sartorij/somas026/CacheSim/traces/429.mcf-184B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/429.mcf-184B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/437.leslie3d-232B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/437.leslie3d-273B.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_74_spp -o MP_mix_74_spp.out -e MP_mix_74_spp.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=spp_dev2 --config=/home/sartorij/somas026/CacheSim/config/spp_dev2.ini  -traces /home/sartorij/somas026/CacheSim/traces/429.mcf-192B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/436.cactusADM-1804B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/481.wrf-455B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/654.roms_s-1613B.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_74_stride -o MP_mix_74_stride.out -e MP_mix_74_stride.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=stride --config=/home/sartorij/somas026/CacheSim/config/stride.ini  -traces /home/sartorij/somas026/CacheSim/traces/429.mcf-192B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/436.cactusADM-1804B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/481.wrf-455B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/654.roms_s-1613B.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_74_bingo -o MP_mix_74_bingo.out -e MP_mix_74_bingo.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=bingo --config=/home/sartorij/somas026/CacheSim/config/bingo.ini  -traces /home/sartorij/somas026/CacheSim/traces/429.mcf-192B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/436.cactusADM-1804B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/481.wrf-455B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/654.roms_s-1613B.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_87_spp -o MP_mix_87_spp.out -e MP_mix_87_spp.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=spp_dev2 --config=/home/sartorij/somas026/CacheSim/config/spp_dev2.ini  -traces /home/sartorij/somas026/CacheSim/traces/429.mcf-51B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/654.roms_s-523B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/602.gcc_s-734B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/459.GemsFDTD-1211B.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_87_mlop -o MP_mix_87_mlop.out -e MP_mix_87_mlop.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=mlop --config=/home/sartorij/somas026/CacheSim/config/mlop.ini  -traces /home/sartorij/somas026/CacheSim/traces/429.mcf-51B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/654.roms_s-523B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/602.gcc_s-734B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/459.GemsFDTD-1211B.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_87_pythia -o MP_mix_87_pythia.out -e MP_mix_87_pythia.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=scooby --config=/home/sartorij/somas026/CacheSim/config/pythia.ini  -traces /home/sartorij/somas026/CacheSim/traces/429.mcf-51B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/654.roms_s-523B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/602.gcc_s-734B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/459.GemsFDTD-1211B.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_92_spp -o MP_mix_92_spp.out -e MP_mix_92_spp.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=spp_dev2 --config=/home/sartorij/somas026/CacheSim/config/spp_dev2.ini  -traces /home/sartorij/somas026/CacheSim/traces/parsec_2.1.streamcluster.simlarge.prebuilt.drop_4750M.length_250M.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/623.xalancbmk_s-202B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/429.mcf-184B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/ligra_Components-Shortcut.com-lj.ungraph.gcc_6.3.0_O3.drop_21000M.length_250M.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_92_stride -o MP_mix_92_stride.out -e MP_mix_92_stride.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=stride --config=/home/sartorij/somas026/CacheSim/config/stride.ini  -traces /home/sartorij/somas026/CacheSim/traces/parsec_2.1.streamcluster.simlarge.prebuilt.drop_4750M.length_250M.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/623.xalancbmk_s-202B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/429.mcf-184B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/ligra_Components-Shortcut.com-lj.ungraph.gcc_6.3.0_O3.drop_21000M.length_250M.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_92_bingo -o MP_mix_92_bingo.out -e MP_mix_92_bingo.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=bingo --config=/home/sartorij/somas026/CacheSim/config/bingo.ini  -traces /home/sartorij/somas026/CacheSim/traces/parsec_2.1.streamcluster.simlarge.prebuilt.drop_4750M.length_250M.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/623.xalancbmk_s-202B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/429.mcf-184B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/ligra_Components-Shortcut.com-lj.ungraph.gcc_6.3.0_O3.drop_21000M.length_250M.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_92_mlop -o MP_mix_92_mlop.out -e MP_mix_92_mlop.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=mlop --config=/home/sartorij/somas026/CacheSim/config/mlop.ini  -traces /home/sartorij/somas026/CacheSim/traces/parsec_2.1.streamcluster.simlarge.prebuilt.drop_4750M.length_250M.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/623.xalancbmk_s-202B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/429.mcf-184B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/ligra_Components-Shortcut.com-lj.ungraph.gcc_6.3.0_O3.drop_21000M.length_250M.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_92_pythia -o MP_mix_92_pythia.out -e MP_mix_92_pythia.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=scooby --config=/home/sartorij/somas026/CacheSim/config/pythia.ini  -traces /home/sartorij/somas026/CacheSim/traces/parsec_2.1.streamcluster.simlarge.prebuilt.drop_4750M.length_250M.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/623.xalancbmk_s-202B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/429.mcf-184B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/ligra_Components-Shortcut.com-lj.ungraph.gcc_6.3.0_O3.drop_21000M.length_250M.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_95_spp -o MP_mix_95_spp.out -e MP_mix_95_spp.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=spp_dev2 --config=/home/sartorij/somas026/CacheSim/config/spp_dev2.ini  -traces /home/sartorij/somas026/CacheSim/traces/ligra_Radii.com-lj.ungraph.gcc_6.3.0_O3.drop_32000M.length_250M.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/ligra_BellmanFord.com-lj.ungraph.gcc_6.3.0_O3.drop_24750M.length_250M.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/654.roms_s-1007B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/ligra_PageRankDelta.com-lj.ungraph.gcc_6.3.0_O3.drop_24000M.length_250M.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_95_mlop -o MP_mix_95_mlop.out -e MP_mix_95_mlop.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=mlop --config=/home/sartorij/somas026/CacheSim/config/mlop.ini  -traces /home/sartorij/somas026/CacheSim/traces/ligra_Radii.com-lj.ungraph.gcc_6.3.0_O3.drop_32000M.length_250M.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/ligra_BellmanFord.com-lj.ungraph.gcc_6.3.0_O3.drop_24750M.length_250M.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/654.roms_s-1007B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/ligra_PageRankDelta.com-lj.ungraph.gcc_6.3.0_O3.drop_24000M.length_250M.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_95_pythia -o MP_mix_95_pythia.out -e MP_mix_95_pythia.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=scooby --config=/home/sartorij/somas026/CacheSim/config/pythia.ini  -traces /home/sartorij/somas026/CacheSim/traces/ligra_Radii.com-lj.ungraph.gcc_6.3.0_O3.drop_32000M.length_250M.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/ligra_BellmanFord.com-lj.ungraph.gcc_6.3.0_O3.drop_24750M.length_250M.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/654.roms_s-1007B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/ligra_PageRankDelta.com-lj.ungraph.gcc_6.3.0_O3.drop_24000M.length_250M.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_101_spp -o MP_mix_101_spp.out -e MP_mix_101_spp.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=spp_dev2 --config=/home/sartorij/somas026/CacheSim/config/spp_dev2.ini  -traces /home/sartorij/somas026/CacheSim/traces/ligra_PageRankDelta.com-lj.ungraph.gcc_6.3.0_O3.drop_17000M.length_250M.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/437.leslie3d-232B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/470.lbm-1274B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/ligra_PageRank.com-lj.ungraph.gcc_6.3.0_O3.drop_79500M.length_250M.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_101_bingo -o MP_mix_101_bingo.out -e MP_mix_101_bingo.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=bingo --config=/home/sartorij/somas026/CacheSim/config/bingo.ini  -traces /home/sartorij/somas026/CacheSim/traces/ligra_PageRankDelta.com-lj.ungraph.gcc_6.3.0_O3.drop_17000M.length_250M.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/437.leslie3d-232B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/470.lbm-1274B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/ligra_PageRank.com-lj.ungraph.gcc_6.3.0_O3.drop_79500M.length_250M.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_101_mlop -o MP_mix_101_mlop.out -e MP_mix_101_mlop.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=mlop --config=/home/sartorij/somas026/CacheSim/config/mlop.ini  -traces /home/sartorij/somas026/CacheSim/traces/ligra_PageRankDelta.com-lj.ungraph.gcc_6.3.0_O3.drop_17000M.length_250M.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/437.leslie3d-232B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/470.lbm-1274B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/ligra_PageRank.com-lj.ungraph.gcc_6.3.0_O3.drop_79500M.length_250M.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_101_dspatch -o MP_mix_101_dspatch.out -e MP_mix_101_dspatch.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=dspatch --config=/home/sartorij/somas026/CacheSim/config/dspatch.ini  -traces /home/sartorij/somas026/CacheSim/traces/ligra_PageRankDelta.com-lj.ungraph.gcc_6.3.0_O3.drop_17000M.length_250M.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/437.leslie3d-232B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/470.lbm-1274B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/ligra_PageRank.com-lj.ungraph.gcc_6.3.0_O3.drop_79500M.length_250M.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_101_pythia -o MP_mix_101_pythia.out -e MP_mix_101_pythia.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=scooby --config=/home/sartorij/somas026/CacheSim/config/pythia.ini  -traces /home/sartorij/somas026/CacheSim/traces/ligra_PageRankDelta.com-lj.ungraph.gcc_6.3.0_O3.drop_17000M.length_250M.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/437.leslie3d-232B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/470.lbm-1274B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/ligra_PageRank.com-lj.ungraph.gcc_6.3.0_O3.drop_79500M.length_250M.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_103_nopref -o MP_mix_103_nopref.out -e MP_mix_103_nopref.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --config=/home/sartorij/somas026/CacheSim/config/nopref.ini  -traces /home/sartorij/somas026/CacheSim/traces/ligra_BellmanFord.com-lj.ungraph.gcc_6.3.0_O3.drop_4000M.length_250M.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/602.gcc_s-2226B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/607.cactuBSSN_s-4004B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/ligra_PageRank.com-lj.ungraph.gcc_6.3.0_O3.drop_60750M.length_250M.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_103_stride -o MP_mix_103_stride.out -e MP_mix_103_stride.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=stride --config=/home/sartorij/somas026/CacheSim/config/stride.ini  -traces /home/sartorij/somas026/CacheSim/traces/ligra_BellmanFord.com-lj.ungraph.gcc_6.3.0_O3.drop_4000M.length_250M.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/602.gcc_s-2226B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/607.cactuBSSN_s-4004B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/ligra_PageRank.com-lj.ungraph.gcc_6.3.0_O3.drop_60750M.length_250M.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_106_spp -o MP_mix_106_spp.out -e MP_mix_106_spp.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=spp_dev2 --config=/home/sartorij/somas026/CacheSim/config/spp_dev2.ini  -traces /home/sartorij/somas026/CacheSim/traces/429.mcf-51B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/ligra_PageRank.com-lj.ungraph.gcc_6.3.0_O3.drop_51000M.length_250M.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/ligra_BFS.com-lj.ungraph.gcc_6.3.0_O3.drop_500M.length_250M.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/603.bwaves_s-1740B.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_106_stride -o MP_mix_106_stride.out -e MP_mix_106_stride.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=stride --config=/home/sartorij/somas026/CacheSim/config/stride.ini  -traces /home/sartorij/somas026/CacheSim/traces/429.mcf-51B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/ligra_PageRank.com-lj.ungraph.gcc_6.3.0_O3.drop_51000M.length_250M.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/ligra_BFS.com-lj.ungraph.gcc_6.3.0_O3.drop_500M.length_250M.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/603.bwaves_s-1740B.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_106_bingo -o MP_mix_106_bingo.out -e MP_mix_106_bingo.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=bingo --config=/home/sartorij/somas026/CacheSim/config/bingo.ini  -traces /home/sartorij/somas026/CacheSim/traces/429.mcf-51B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/ligra_PageRank.com-lj.ungraph.gcc_6.3.0_O3.drop_51000M.length_250M.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/ligra_BFS.com-lj.ungraph.gcc_6.3.0_O3.drop_500M.length_250M.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/603.bwaves_s-1740B.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_106_mlop -o MP_mix_106_mlop.out -e MP_mix_106_mlop.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=mlop --config=/home/sartorij/somas026/CacheSim/config/mlop.ini  -traces /home/sartorij/somas026/CacheSim/traces/429.mcf-51B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/ligra_PageRank.com-lj.ungraph.gcc_6.3.0_O3.drop_51000M.length_250M.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/ligra_BFS.com-lj.ungraph.gcc_6.3.0_O3.drop_500M.length_250M.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/603.bwaves_s-1740B.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_106_pythia -o MP_mix_106_pythia.out -e MP_mix_106_pythia.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=scooby --config=/home/sartorij/somas026/CacheSim/config/pythia.ini  -traces /home/sartorij/somas026/CacheSim/traces/429.mcf-51B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/ligra_PageRank.com-lj.ungraph.gcc_6.3.0_O3.drop_51000M.length_250M.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/ligra_BFS.com-lj.ungraph.gcc_6.3.0_O3.drop_500M.length_250M.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/603.bwaves_s-1740B.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_126_spp -o MP_mix_126_spp.out -e MP_mix_126_spp.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=spp_dev2 --config=/home/sartorij/somas026/CacheSim/config/spp_dev2.ini  -traces /home/sartorij/somas026/CacheSim/traces/ligra_BFSCC.com-lj.ungraph.gcc_6.3.0_O3.drop_18750M.length_250M.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/ligra_PageRankDelta.com-lj.ungraph.gcc_6.3.0_O3.drop_52000M.length_250M.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/437.leslie3d-273B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/429.mcf-51B.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_126_bingo -o MP_mix_126_bingo.out -e MP_mix_126_bingo.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=bingo --config=/home/sartorij/somas026/CacheSim/config/bingo.ini  -traces /home/sartorij/somas026/CacheSim/traces/ligra_BFSCC.com-lj.ungraph.gcc_6.3.0_O3.drop_18750M.length_250M.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/ligra_PageRankDelta.com-lj.ungraph.gcc_6.3.0_O3.drop_52000M.length_250M.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/437.leslie3d-273B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/429.mcf-51B.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_128_mlop -o MP_mix_128_mlop.out -e MP_mix_128_mlop.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=mlop --config=/home/sartorij/somas026/CacheSim/config/mlop.ini  -traces /home/sartorij/somas026/CacheSim/traces/parsec_2.1.canneal.simlarge.prebuilt.drop_3000M.length_250M.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/ligra_Triangle.com-lj.ungraph.gcc_6.3.0_O3.drop_18000M.length_250M.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/619.lbm_s-2677B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/620.omnetpp_s-874B.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_128_dspatch -o MP_mix_128_dspatch.out -e MP_mix_128_dspatch.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=dspatch --config=/home/sartorij/somas026/CacheSim/config/dspatch.ini  -traces /home/sartorij/somas026/CacheSim/traces/parsec_2.1.canneal.simlarge.prebuilt.drop_3000M.length_250M.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/ligra_Triangle.com-lj.ungraph.gcc_6.3.0_O3.drop_18000M.length_250M.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/619.lbm_s-2677B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/620.omnetpp_s-874B.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_134_nopref -o MP_mix_134_nopref.out -e MP_mix_134_nopref.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --config=/home/sartorij/somas026/CacheSim/config/nopref.ini  -traces /home/sartorij/somas026/CacheSim/traces/470.lbm-1274B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/429.mcf-192B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/ligra_BFS.com-lj.ungraph.gcc_6.3.0_O3.drop_21500M.length_250M.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/462.libquantum-714B.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_134_spp -o MP_mix_134_spp.out -e MP_mix_134_spp.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=spp_dev2 --config=/home/sartorij/somas026/CacheSim/config/spp_dev2.ini  -traces /home/sartorij/somas026/CacheSim/traces/470.lbm-1274B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/429.mcf-192B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/ligra_BFS.com-lj.ungraph.gcc_6.3.0_O3.drop_21500M.length_250M.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/462.libquantum-714B.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_134_stride -o MP_mix_134_stride.out -e MP_mix_134_stride.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=stride --config=/home/sartorij/somas026/CacheSim/config/stride.ini  -traces /home/sartorij/somas026/CacheSim/traces/470.lbm-1274B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/429.mcf-192B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/ligra_BFS.com-lj.ungraph.gcc_6.3.0_O3.drop_21500M.length_250M.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/462.libquantum-714B.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_134_bingo -o MP_mix_134_bingo.out -e MP_mix_134_bingo.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=bingo --config=/home/sartorij/somas026/CacheSim/config/bingo.ini  -traces /home/sartorij/somas026/CacheSim/traces/470.lbm-1274B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/429.mcf-192B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/ligra_BFS.com-lj.ungraph.gcc_6.3.0_O3.drop_21500M.length_250M.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/462.libquantum-714B.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_134_mlop -o MP_mix_134_mlop.out -e MP_mix_134_mlop.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=mlop --config=/home/sartorij/somas026/CacheSim/config/mlop.ini  -traces /home/sartorij/somas026/CacheSim/traces/470.lbm-1274B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/429.mcf-192B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/ligra_BFS.com-lj.ungraph.gcc_6.3.0_O3.drop_21500M.length_250M.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/462.libquantum-714B.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_134_dspatch -o MP_mix_134_dspatch.out -e MP_mix_134_dspatch.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=dspatch --config=/home/sartorij/somas026/CacheSim/config/dspatch.ini  -traces /home/sartorij/somas026/CacheSim/traces/470.lbm-1274B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/429.mcf-192B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/ligra_BFS.com-lj.ungraph.gcc_6.3.0_O3.drop_21500M.length_250M.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/462.libquantum-714B.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_134_pythia -o MP_mix_134_pythia.out -e MP_mix_134_pythia.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=scooby --config=/home/sartorij/somas026/CacheSim/config/pythia.ini  -traces /home/sartorij/somas026/CacheSim/traces/470.lbm-1274B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/429.mcf-192B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/ligra_BFS.com-lj.ungraph.gcc_6.3.0_O3.drop_21500M.length_250M.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/462.libquantum-714B.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_152_mlop -o MP_mix_152_mlop.out -e MP_mix_152_mlop.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=mlop --config=/home/sartorij/somas026/CacheSim/config/mlop.ini  -traces /home/sartorij/somas026/CacheSim/traces/parsec_2.1.canneal.simlarge.prebuilt.drop_4750M.length_250M.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/ligra_Components-Shortcut.com-lj.ungraph.gcc_6.3.0_O3.drop_3500M.length_250M.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/445.gobmk-36B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/619.lbm_s-2677B.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_179_pythia -o MP_mix_179_pythia.out -e MP_mix_179_pythia.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=scooby --config=/home/sartorij/somas026/CacheSim/config/pythia.ini  -traces /home/sartorij/somas026/CacheSim/traces/605.mcf_s-1536B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/607.cactuBSSN_s-2421B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/ligra_PageRank.com-lj.ungraph.gcc_6.3.0_O3.drop_51000M.length_250M.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/ligra_Components.com-lj.ungraph.gcc_6.3.0_O3.drop_23750M.length_250M.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_183_spp -o MP_mix_183_spp.out -e MP_mix_183_spp.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=spp_dev2 --config=/home/sartorij/somas026/CacheSim/config/spp_dev2.ini  -traces /home/sartorij/somas026/CacheSim/traces/ligra_Components-Shortcut.com-lj.ungraph.gcc_6.3.0_O3.drop_750M.length_250M.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/429.mcf-51B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/ligra_BFSCC.com-lj.ungraph.gcc_6.3.0_O3.drop_3500M.length_250M.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/ligra_MIS.com-lj.ungraph.gcc_6.3.0_O3.drop_5000M.length_250M.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_183_stride -o MP_mix_183_stride.out -e MP_mix_183_stride.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=stride --config=/home/sartorij/somas026/CacheSim/config/stride.ini  -traces /home/sartorij/somas026/CacheSim/traces/ligra_Components-Shortcut.com-lj.ungraph.gcc_6.3.0_O3.drop_750M.length_250M.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/429.mcf-51B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/ligra_BFSCC.com-lj.ungraph.gcc_6.3.0_O3.drop_3500M.length_250M.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/ligra_MIS.com-lj.ungraph.gcc_6.3.0_O3.drop_5000M.length_250M.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_183_bingo -o MP_mix_183_bingo.out -e MP_mix_183_bingo.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=bingo --config=/home/sartorij/somas026/CacheSim/config/bingo.ini  -traces /home/sartorij/somas026/CacheSim/traces/ligra_Components-Shortcut.com-lj.ungraph.gcc_6.3.0_O3.drop_750M.length_250M.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/429.mcf-51B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/ligra_BFSCC.com-lj.ungraph.gcc_6.3.0_O3.drop_3500M.length_250M.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/ligra_MIS.com-lj.ungraph.gcc_6.3.0_O3.drop_5000M.length_250M.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_183_mlop -o MP_mix_183_mlop.out -e MP_mix_183_mlop.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=mlop --config=/home/sartorij/somas026/CacheSim/config/mlop.ini  -traces /home/sartorij/somas026/CacheSim/traces/ligra_Components-Shortcut.com-lj.ungraph.gcc_6.3.0_O3.drop_750M.length_250M.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/429.mcf-51B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/ligra_BFSCC.com-lj.ungraph.gcc_6.3.0_O3.drop_3500M.length_250M.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/ligra_MIS.com-lj.ungraph.gcc_6.3.0_O3.drop_5000M.length_250M.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_188_bingo -o MP_mix_188_bingo.out -e MP_mix_188_bingo.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=bingo --config=/home/sartorij/somas026/CacheSim/config/bingo.ini  -traces /home/sartorij/somas026/CacheSim/traces/481.wrf-455B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/ligra_PageRankDelta.com-lj.ungraph.gcc_6.3.0_O3.drop_24000M.length_250M.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/619.lbm_s-2676B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/ligra_PageRank.com-lj.ungraph.gcc_6.3.0_O3.drop_21750M.length_250M.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_188_mlop -o MP_mix_188_mlop.out -e MP_mix_188_mlop.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=mlop --config=/home/sartorij/somas026/CacheSim/config/mlop.ini  -traces /home/sartorij/somas026/CacheSim/traces/481.wrf-455B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/ligra_PageRankDelta.com-lj.ungraph.gcc_6.3.0_O3.drop_24000M.length_250M.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/619.lbm_s-2676B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/ligra_PageRank.com-lj.ungraph.gcc_6.3.0_O3.drop_21750M.length_250M.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_199_nopref -o MP_mix_199_nopref.out -e MP_mix_199_nopref.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --config=/home/sartorij/somas026/CacheSim/config/nopref.ini  -traces /home/sartorij/somas026/CacheSim/traces/ligra_BFS-Bitvector.com-lj.ungraph.gcc_6.3.0_O3.drop_500M.length_250M.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/481.wrf-816B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/ligra_Components.com-lj.ungraph.gcc_6.3.0_O3.drop_23750M.length_250M.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/429.mcf-192B.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_199_spp -o MP_mix_199_spp.out -e MP_mix_199_spp.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=spp_dev2 --config=/home/sartorij/somas026/CacheSim/config/spp_dev2.ini  -traces /home/sartorij/somas026/CacheSim/traces/ligra_BFS-Bitvector.com-lj.ungraph.gcc_6.3.0_O3.drop_500M.length_250M.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/481.wrf-816B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/ligra_Components.com-lj.ungraph.gcc_6.3.0_O3.drop_23750M.length_250M.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/429.mcf-192B.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_199_stride -o MP_mix_199_stride.out -e MP_mix_199_stride.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=stride --config=/home/sartorij/somas026/CacheSim/config/stride.ini  -traces /home/sartorij/somas026/CacheSim/traces/ligra_BFS-Bitvector.com-lj.ungraph.gcc_6.3.0_O3.drop_500M.length_250M.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/481.wrf-816B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/ligra_Components.com-lj.ungraph.gcc_6.3.0_O3.drop_23750M.length_250M.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/429.mcf-192B.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_199_bingo -o MP_mix_199_bingo.out -e MP_mix_199_bingo.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=bingo --config=/home/sartorij/somas026/CacheSim/config/bingo.ini  -traces /home/sartorij/somas026/CacheSim/traces/ligra_BFS-Bitvector.com-lj.ungraph.gcc_6.3.0_O3.drop_500M.length_250M.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/481.wrf-816B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/ligra_Components.com-lj.ungraph.gcc_6.3.0_O3.drop_23750M.length_250M.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/429.mcf-192B.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_199_mlop -o MP_mix_199_mlop.out -e MP_mix_199_mlop.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=mlop --config=/home/sartorij/somas026/CacheSim/config/mlop.ini  -traces /home/sartorij/somas026/CacheSim/traces/ligra_BFS-Bitvector.com-lj.ungraph.gcc_6.3.0_O3.drop_500M.length_250M.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/481.wrf-816B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/ligra_Components.com-lj.ungraph.gcc_6.3.0_O3.drop_23750M.length_250M.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/429.mcf-192B.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_199_pythia -o MP_mix_199_pythia.out -e MP_mix_199_pythia.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=scooby --config=/home/sartorij/somas026/CacheSim/config/pythia.ini  -traces /home/sartorij/somas026/CacheSim/traces/ligra_BFS-Bitvector.com-lj.ungraph.gcc_6.3.0_O3.drop_500M.length_250M.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/481.wrf-816B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/ligra_Components.com-lj.ungraph.gcc_6.3.0_O3.drop_23750M.length_250M.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/429.mcf-192B.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_244_spp -o MP_mix_244_spp.out -e MP_mix_244_spp.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=spp_dev2 --config=/home/sartorij/somas026/CacheSim/config/spp_dev2.ini  -traces /home/sartorij/somas026/CacheSim/traces/623.xalancbmk_s-10B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/ligra_BFSCC.com-lj.ungraph.gcc_6.3.0_O3.drop_15500M.length_250M.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/ligra_Radii.com-lj.ungraph.gcc_6.3.0_O3.drop_18000M.length_250M.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/ligra_MIS.com-lj.ungraph.gcc_6.3.0_O3.drop_3500M.length_250M.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_244_bingo -o MP_mix_244_bingo.out -e MP_mix_244_bingo.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=bingo --config=/home/sartorij/somas026/CacheSim/config/bingo.ini  -traces /home/sartorij/somas026/CacheSim/traces/623.xalancbmk_s-10B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/ligra_BFSCC.com-lj.ungraph.gcc_6.3.0_O3.drop_15500M.length_250M.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/ligra_Radii.com-lj.ungraph.gcc_6.3.0_O3.drop_18000M.length_250M.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/ligra_MIS.com-lj.ungraph.gcc_6.3.0_O3.drop_3500M.length_250M.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_244_mlop -o MP_mix_244_mlop.out -e MP_mix_244_mlop.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=mlop --config=/home/sartorij/somas026/CacheSim/config/mlop.ini  -traces /home/sartorij/somas026/CacheSim/traces/623.xalancbmk_s-10B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/ligra_BFSCC.com-lj.ungraph.gcc_6.3.0_O3.drop_15500M.length_250M.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/ligra_Radii.com-lj.ungraph.gcc_6.3.0_O3.drop_18000M.length_250M.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/ligra_MIS.com-lj.ungraph.gcc_6.3.0_O3.drop_3500M.length_250M.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_255_mlop -o MP_mix_255_mlop.out -e MP_mix_255_mlop.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=mlop --config=/home/sartorij/somas026/CacheSim/config/mlop.ini  -traces /home/sartorij/somas026/CacheSim/traces/620.omnetpp_s-141B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/654.roms_s-523B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/603.bwaves_s-2609B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/parsec_2.1.canneal.simlarge.prebuilt.drop_4750M.length_250M.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_258_nopref -o MP_mix_258_nopref.out -e MP_mix_258_nopref.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --config=/home/sartorij/somas026/CacheSim/config/nopref.ini  -traces /home/sartorij/somas026/CacheSim/traces/429.mcf-192B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/ligra_BC.com-lj.ungraph.gcc_6.3.0_O3.drop_3500M.length_250M.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/456.hmmer-327B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/ligra_CF.com-lj.ungraph.gcc_6.3.0_O3.drop_2500M.length_250M.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_258_spp -o MP_mix_258_spp.out -e MP_mix_258_spp.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=spp_dev2 --config=/home/sartorij/somas026/CacheSim/config/spp_dev2.ini  -traces /home/sartorij/somas026/CacheSim/traces/429.mcf-192B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/ligra_BC.com-lj.ungraph.gcc_6.3.0_O3.drop_3500M.length_250M.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/456.hmmer-327B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/ligra_CF.com-lj.ungraph.gcc_6.3.0_O3.drop_2500M.length_250M.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_258_stride -o MP_mix_258_stride.out -e MP_mix_258_stride.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=stride --config=/home/sartorij/somas026/CacheSim/config/stride.ini  -traces /home/sartorij/somas026/CacheSim/traces/429.mcf-192B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/ligra_BC.com-lj.ungraph.gcc_6.3.0_O3.drop_3500M.length_250M.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/456.hmmer-327B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/ligra_CF.com-lj.ungraph.gcc_6.3.0_O3.drop_2500M.length_250M.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_258_bingo -o MP_mix_258_bingo.out -e MP_mix_258_bingo.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=bingo --config=/home/sartorij/somas026/CacheSim/config/bingo.ini  -traces /home/sartorij/somas026/CacheSim/traces/429.mcf-192B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/ligra_BC.com-lj.ungraph.gcc_6.3.0_O3.drop_3500M.length_250M.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/456.hmmer-327B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/ligra_CF.com-lj.ungraph.gcc_6.3.0_O3.drop_2500M.length_250M.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_258_mlop -o MP_mix_258_mlop.out -e MP_mix_258_mlop.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=mlop --config=/home/sartorij/somas026/CacheSim/config/mlop.ini  -traces /home/sartorij/somas026/CacheSim/traces/429.mcf-192B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/ligra_BC.com-lj.ungraph.gcc_6.3.0_O3.drop_3500M.length_250M.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/456.hmmer-327B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/ligra_CF.com-lj.ungraph.gcc_6.3.0_O3.drop_2500M.length_250M.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_274_spp -o MP_mix_274_spp.out -e MP_mix_274_spp.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=spp_dev2 --config=/home/sartorij/somas026/CacheSim/config/spp_dev2.ini  -traces /home/sartorij/somas026/CacheSim/traces/ligra_CF.com-lj.ungraph.gcc_6.3.0_O3.drop_184750M.length_250M.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/ligra_BFS-Bitvector.com-lj.ungraph.gcc_6.3.0_O3.drop_18000M.length_250M.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/605.mcf_s-782B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/ligra_BFSCC.com-lj.ungraph.gcc_6.3.0_O3.drop_17000M.length_250M.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_274_mlop -o MP_mix_274_mlop.out -e MP_mix_274_mlop.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=mlop --config=/home/sartorij/somas026/CacheSim/config/mlop.ini  -traces /home/sartorij/somas026/CacheSim/traces/ligra_CF.com-lj.ungraph.gcc_6.3.0_O3.drop_184750M.length_250M.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/ligra_BFS-Bitvector.com-lj.ungraph.gcc_6.3.0_O3.drop_18000M.length_250M.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/605.mcf_s-782B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/ligra_BFSCC.com-lj.ungraph.gcc_6.3.0_O3.drop_17000M.length_250M.champsimtrace.xz"
sbatch -p agsmall --mincpus=1 -t 5:00:00 -c 4 -J MP_mix_274_pythia -o MP_mix_274_pythia.out -e MP_mix_274_pythia.err /home/sartorij/somas026/CacheSim/wrapper.sh /home/sartorij/somas026/CacheSim/bin/perceptron-multi-multi-no-srrip-4core "--warmup_instructions=50000000 --simulation_instructions=150000000 --l2c_prefetcher_types=scooby --config=/home/sartorij/somas026/CacheSim/config/pythia.ini  -traces /home/sartorij/somas026/CacheSim/traces/ligra_CF.com-lj.ungraph.gcc_6.3.0_O3.drop_184750M.length_250M.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/ligra_BFS-Bitvector.com-lj.ungraph.gcc_6.3.0_O3.drop_18000M.length_250M.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/605.mcf_s-782B.champsimtrace.xz /home/sartorij/somas026/CacheSim/traces/ligra_BFSCC.com-lj.ungraph.gcc_6.3.0_O3.drop_17000M.length_250M.champsimtrace.xz"
