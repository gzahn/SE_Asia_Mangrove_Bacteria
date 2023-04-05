#!/bin/bash
set -euxo pipefail

if ! command -v prefetch &> /dev/null
then
  echo "Please install the sra-toolkit"
  echo "https://hpc.nih.gov/apps/sratoolkit.html"
  exit 1
else
  echo "sra-toolkit found in PATH. Proceeding to download raw data from PRJNA735404..."
fi

# download SRA files
prefetch PRJNA735404



# extract raw fastq files
fasterq-dump SRR14740814
fasterq-dump SRR14740815
fasterq-dump SRR14740816
fasterq-dump SRR14740817
fasterq-dump SRR14740818
fasterq-dump SRR14740819
fasterq-dump SRR14740820
fasterq-dump SRR14740821
fasterq-dump SRR14740822
fasterq-dump SRR14740823
fasterq-dump SRR14740824
fasterq-dump SRR14740825
fasterq-dump SRR14740826
fasterq-dump SRR14740827
fasterq-dump SRR14740828
fasterq-dump SRR14740829
fasterq-dump SRR14740830
fasterq-dump SRR14740831
fasterq-dump SRR14740832
fasterq-dump SRR14740833
fasterq-dump SRR14740834
fasterq-dump SRR14740835
fasterq-dump SRR14740836
fasterq-dump SRR14740837
fasterq-dump SRR14740838
fasterq-dump SRR14740839
fasterq-dump SRR14740840
fasterq-dump SRR14740841
fasterq-dump SRR14740842
fasterq-dump SRR14740843
fasterq-dump SRR14740844
fasterq-dump SRR14740845
fasterq-dump SRR14740846
fasterq-dump SRR14740847
fasterq-dump SRR14740848
fasterq-dump SRR14740849
fasterq-dump SRR14740850
fasterq-dump SRR14740851
fasterq-dump SRR14740852
fasterq-dump SRR14740853
fasterq-dump SRR14740854
fasterq-dump SRR14740855
fasterq-dump SRR14740856
fasterq-dump SRR14740857
fasterq-dump SRR14740858
fasterq-dump SRR14740859
fasterq-dump SRR14740860
fasterq-dump SRR14740861
fasterq-dump SRR14740862
fasterq-dump SRR14740863
fasterq-dump SRR14740864
fasterq-dump SRR14740865
fasterq-dump SRR14740866
fasterq-dump SRR14740867
fasterq-dump SRR14740868
fasterq-dump SRR14740869
fasterq-dump SRR14740870
fasterq-dump SRR14740871
fasterq-dump SRR14740872
fasterq-dump SRR14740873
fasterq-dump SRR14740874
fasterq-dump SRR14740875
fasterq-dump SRR14740876
fasterq-dump SRR14740877
fasterq-dump SRR14740878
fasterq-dump SRR14740879
fasterq-dump SRR14740880
fasterq-dump SRR14740881
fasterq-dump SRR14740882
fasterq-dump SRR14740883
fasterq-dump SRR14740884
fasterq-dump SRR14740885
fasterq-dump SRR14740886
fasterq-dump SRR14740887
fasterq-dump SRR14740888
fasterq-dump SRR14740889
fasterq-dump SRR14740890
fasterq-dump SRR14740891
fasterq-dump SRR14740892
fasterq-dump SRR14740893
fasterq-dump SRR14740894
fasterq-dump SRR14740895
fasterq-dump SRR14740896
fasterq-dump SRR14740897
fasterq-dump SRR14740898
fasterq-dump SRR14740899
fasterq-dump SRR14740900
fasterq-dump SRR14740901
fasterq-dump SRR14740902
fasterq-dump SRR14740903
fasterq-dump SRR14740904
fasterq-dump SRR14740905
fasterq-dump SRR14740906
fasterq-dump SRR14740907
fasterq-dump SRR14740908
fasterq-dump SRR14740909
fasterq-dump SRR14740910
fasterq-dump SRR14740911
fasterq-dump SRR14740912
fasterq-dump SRR14740913
fasterq-dump SRR14740914
fasterq-dump SRR14740915
fasterq-dump SRR14740916
fasterq-dump SRR14740917
fasterq-dump SRR14740918
fasterq-dump SRR14740919
fasterq-dump SRR14740920
fasterq-dump SRR14740921
fasterq-dump SRR14740922
fasterq-dump SRR14740923
fasterq-dump SRR14740924
fasterq-dump SRR14740925
fasterq-dump SRR14740926
fasterq-dump SRR14740927
fasterq-dump SRR14740928
fasterq-dump SRR14740929
fasterq-dump SRR14740930
fasterq-dump SRR14740931
fasterq-dump SRR14740932
fasterq-dump SRR14740933
fasterq-dump SRR14740934
fasterq-dump SRR14740935
fasterq-dump SRR14740936
fasterq-dump SRR14740937
fasterq-dump SRR14740938
fasterq-dump SRR14740939
fasterq-dump SRR14740940
fasterq-dump SRR14740941
fasterq-dump SRR14740942
fasterq-dump SRR14740943
fasterq-dump SRR14740944
fasterq-dump SRR14740945
fasterq-dump SRR14740946
fasterq-dump SRR14740947
fasterq-dump SRR14740948
fasterq-dump SRR14740949
fasterq-dump SRR14740950
fasterq-dump SRR14740951
fasterq-dump SRR14740952
fasterq-dump SRR14740953
fasterq-dump SRR14740954
fasterq-dump SRR14740955
fasterq-dump SRR14740956
fasterq-dump SRR14740957
fasterq-dump SRR14740958
fasterq-dump SRR14740959
fasterq-dump SRR14740960
fasterq-dump SRR14740961
fasterq-dump SRR14740962
fasterq-dump SRR14740963
fasterq-dump SRR14740964
fasterq-dump SRR14740965
fasterq-dump SRR14740966
fasterq-dump SRR14740967
fasterq-dump SRR14740968
fasterq-dump SRR14740969
fasterq-dump SRR14740970
fasterq-dump SRR14740971
fasterq-dump SRR14740972
fasterq-dump SRR14740973
fasterq-dump SRR14740974
fasterq-dump SRR14740975
fasterq-dump SRR14740976
fasterq-dump SRR14740977
fasterq-dump SRR14740978
fasterq-dump SRR14740979
fasterq-dump SRR14740980
fasterq-dump SRR14740981
fasterq-dump SRR14740982
fasterq-dump SRR14740983
fasterq-dump SRR14740984
fasterq-dump SRR14740985
fasterq-dump SRR14740986
fasterq-dump SRR14740987
fasterq-dump SRR14740988
fasterq-dump SRR14740989
fasterq-dump SRR14740990
fasterq-dump SRR14740991
fasterq-dump SRR14740992
fasterq-dump SRR14740993
fasterq-dump SRR14740994
fasterq-dump SRR14740995
fasterq-dump SRR14740996
fasterq-dump SRR14740997
fasterq-dump SRR14740998
fasterq-dump SRR14740999
fasterq-dump SRR14741000
fasterq-dump SRR14741001
fasterq-dump SRR14741002
fasterq-dump SRR14741003
fasterq-dump SRR14741004
fasterq-dump SRR14741005
fasterq-dump SRR14741006
fasterq-dump SRR14741007
fasterq-dump SRR14741008
fasterq-dump SRR14741009
fasterq-dump SRR14741010
fasterq-dump SRR14741011
fasterq-dump SRR14741012
fasterq-dump SRR14741013
fasterq-dump SRR14741014
fasterq-dump SRR14741015
fasterq-dump SRR14741016
fasterq-dump SRR14741017
fasterq-dump SRR14741018
fasterq-dump SRR14741019
fasterq-dump SRR14741020
fasterq-dump SRR14741021
fasterq-dump SRR14741022
fasterq-dump SRR14741023
fasterq-dump SRR14741024
fasterq-dump SRR14741025
fasterq-dump SRR14741026
fasterq-dump SRR14741027
fasterq-dump SRR14741028
fasterq-dump SRR14741029
fasterq-dump SRR14741030
fasterq-dump SRR14741031
fasterq-dump SRR14741032
fasterq-dump SRR14741033
fasterq-dump SRR14741034
fasterq-dump SRR14741035
fasterq-dump SRR14741036
fasterq-dump SRR14741037
fasterq-dump SRR14741038
fasterq-dump SRR14741039
fasterq-dump SRR14741040
fasterq-dump SRR14741041
fasterq-dump SRR14741042
fasterq-dump SRR14741043
fasterq-dump SRR14741044
fasterq-dump SRR14741045
fasterq-dump SRR14741046
fasterq-dump SRR14741047
fasterq-dump SRR14741048
fasterq-dump SRR14741049
fasterq-dump SRR14741050
fasterq-dump SRR14741051
fasterq-dump SRR14741052
fasterq-dump SRR14741053
fasterq-dump SRR14741054
fasterq-dump SRR14741055
fasterq-dump SRR14741056
fasterq-dump SRR14741057
fasterq-dump SRR14741058
fasterq-dump SRR14741059
fasterq-dump SRR14741060
fasterq-dump SRR14741061
fasterq-dump SRR14741062
fasterq-dump SRR14741063
fasterq-dump SRR14741064
fasterq-dump SRR14741065
fasterq-dump SRR14741066
fasterq-dump SRR14741067
fasterq-dump SRR14741068
fasterq-dump SRR14741069
fasterq-dump SRR14741070
fasterq-dump SRR14741071
fasterq-dump SRR14741072
fasterq-dump SRR14741073
fasterq-dump SRR14741074
fasterq-dump SRR14741075
fasterq-dump SRR14741076
fasterq-dump SRR14741077
fasterq-dump SRR14741078
fasterq-dump SRR14741079
fasterq-dump SRR14741080
fasterq-dump SRR14741081
fasterq-dump SRR14741082
fasterq-dump SRR14741083
fasterq-dump SRR14741084
fasterq-dump SRR14741085
fasterq-dump SRR14741086
fasterq-dump SRR14741087
fasterq-dump SRR14741088
fasterq-dump SRR14741089
fasterq-dump SRR14741090
fasterq-dump SRR14741091
fasterq-dump SRR14741092
fasterq-dump SRR14741093
fasterq-dump SRR14741094
fasterq-dump SRR14741095
fasterq-dump SRR14741096
fasterq-dump SRR14741097
fasterq-dump SRR14741098
fasterq-dump SRR14741099
fasterq-dump SRR14741100
fasterq-dump SRR14741101
fasterq-dump SRR14741102
fasterq-dump SRR14741103
fasterq-dump SRR14741104
fasterq-dump SRR14741105
fasterq-dump SRR14741106
fasterq-dump SRR14741107
fasterq-dump SRR14741108
fasterq-dump SRR14741109
fasterq-dump SRR14741110
fasterq-dump SRR14741111
fasterq-dump SRR14741112
fasterq-dump SRR14741113
fasterq-dump SRR14741114
fasterq-dump SRR14741115
fasterq-dump SRR14741116
fasterq-dump SRR14741117
fasterq-dump SRR14741118
fasterq-dump SRR14741119
fasterq-dump SRR14741120
fasterq-dump SRR14741121
fasterq-dump SRR14741122
fasterq-dump SRR14741123
fasterq-dump SRR14741124
fasterq-dump SRR14741125
fasterq-dump SRR14741126
fasterq-dump SRR14741127
fasterq-dump SRR14741128
fasterq-dump SRR14741129
fasterq-dump SRR14741130
fasterq-dump SRR14741131
fasterq-dump SRR14741132
fasterq-dump SRR14741133
fasterq-dump SRR14741134
fasterq-dump SRR14741135
fasterq-dump SRR14741136
fasterq-dump SRR14741137
fasterq-dump SRR14741138
fasterq-dump SRR14741139
fasterq-dump SRR14741140
fasterq-dump SRR14741141
fasterq-dump SRR14741142
fasterq-dump SRR14741143
fasterq-dump SRR14741144
fasterq-dump SRR14741145
fasterq-dump SRR14741146
fasterq-dump SRR14741147
fasterq-dump SRR14741148
fasterq-dump SRR14741149
fasterq-dump SRR14741150
fasterq-dump SRR14741151
fasterq-dump SRR14741152
fasterq-dump SRR14741153
fasterq-dump SRR14741154
fasterq-dump SRR14741155
fasterq-dump SRR14741156
fasterq-dump SRR14741157
fasterq-dump SRR14741158
fasterq-dump SRR14741159
fasterq-dump SRR14741160
fasterq-dump SRR14741161
fasterq-dump SRR14741162
fasterq-dump SRR14741163
fasterq-dump SRR14741164
fasterq-dump SRR14741165
fasterq-dump SRR14741166
fasterq-dump SRR14741167
fasterq-dump SRR14741168
fasterq-dump SRR14741169
fasterq-dump SRR14741170
fasterq-dump SRR14741171
fasterq-dump SRR14741172
fasterq-dump SRR14741173
fasterq-dump SRR14741174
fasterq-dump SRR14741175
fasterq-dump SRR14741176
fasterq-dump SRR14741177
fasterq-dump SRR14741178
fasterq-dump SRR14741179
fasterq-dump SRR14741180
fasterq-dump SRR14741181
fasterq-dump SRR14741182
fasterq-dump SRR14741183
fasterq-dump SRR14741184
fasterq-dump SRR14741185
fasterq-dump SRR14741186
fasterq-dump SRR14741187
fasterq-dump SRR14741188
fasterq-dump SRR14741189
fasterq-dump SRR14741190
fasterq-dump SRR14741191
fasterq-dump SRR14741192
fasterq-dump SRR14741193
fasterq-dump SRR14741194
fasterq-dump SRR14741195
fasterq-dump SRR14741196
fasterq-dump SRR14741197
fasterq-dump SRR14741198
fasterq-dump SRR14741199
fasterq-dump SRR14741200
fasterq-dump SRR14741201
fasterq-dump SRR14741202
fasterq-dump SRR14741203
fasterq-dump SRR14741204
fasterq-dump SRR14741205
fasterq-dump SRR14741206
fasterq-dump SRR14741207
fasterq-dump SRR14741208
fasterq-dump SRR14741209
fasterq-dump SRR14741210
fasterq-dump SRR14741211
fasterq-dump SRR14741212
fasterq-dump SRR14741213
fasterq-dump SRR14741214
fasterq-dump SRR14741215
fasterq-dump SRR14741216
fasterq-dump SRR14741217
fasterq-dump SRR14741218
fasterq-dump SRR14741219
fasterq-dump SRR14741220
fasterq-dump SRR14741221
fasterq-dump SRR14741222
fasterq-dump SRR14741223
fasterq-dump SRR14741224
fasterq-dump SRR14741225
fasterq-dump SRR14741226
fasterq-dump SRR14741227
fasterq-dump SRR14741228
fasterq-dump SRR14741229
fasterq-dump SRR14741230
fasterq-dump SRR14741231
fasterq-dump SRR14741232
fasterq-dump SRR14741233
fasterq-dump SRR14741234
fasterq-dump SRR14741235
fasterq-dump SRR14741236
fasterq-dump SRR14741237
fasterq-dump SRR14741238
fasterq-dump SRR14741239
fasterq-dump SRR14741240
fasterq-dump SRR14741241
fasterq-dump SRR14741242
fasterq-dump SRR14741243
fasterq-dump SRR14741244
fasterq-dump SRR14741245
fasterq-dump SRR14741246
fasterq-dump SRR14741247
fasterq-dump SRR14741248
fasterq-dump SRR14741249
fasterq-dump SRR14741250
fasterq-dump SRR14741251
fasterq-dump SRR14741252
fasterq-dump SRR14741253
fasterq-dump SRR14741254
fasterq-dump SRR14741255
fasterq-dump SRR14741256
fasterq-dump SRR14741257
fasterq-dump SRR14741258
fasterq-dump SRR14741259
fasterq-dump SRR14741260
fasterq-dump SRR14741261
fasterq-dump SRR14741262
fasterq-dump SRR14741263
fasterq-dump SRR14741264
fasterq-dump SRR14741265
fasterq-dump SRR14741266
fasterq-dump SRR14741267
fasterq-dump SRR14741268
fasterq-dump SRR14741269
fasterq-dump SRR14741270
fasterq-dump SRR14741271
fasterq-dump SRR14741272
fasterq-dump SRR14741273
fasterq-dump SRR14741274
fasterq-dump SRR14741275
fasterq-dump SRR14741276
fasterq-dump SRR14741277
fasterq-dump SRR14741278
fasterq-dump SRR14741279
fasterq-dump SRR14741280
fasterq-dump SRR14741281
fasterq-dump SRR14741282
fasterq-dump SRR14741283
fasterq-dump SRR14741284
fasterq-dump SRR14741285
fasterq-dump SRR14741286
fasterq-dump SRR14741287
fasterq-dump SRR14741288
fasterq-dump SRR14741289
fasterq-dump SRR14741290
fasterq-dump SRR14741291
fasterq-dump SRR14741292
fasterq-dump SRR14741293
fasterq-dump SRR14741294
fasterq-dump SRR14741295
fasterq-dump SRR14741296
fasterq-dump SRR14741297
fasterq-dump SRR14741298
fasterq-dump SRR14741299
fasterq-dump SRR14741300
fasterq-dump SRR14741301
fasterq-dump SRR14741302
fasterq-dump SRR14741303
fasterq-dump SRR14741304
fasterq-dump SRR14741305
fasterq-dump SRR14741306
fasterq-dump SRR14741307
fasterq-dump SRR14741308
fasterq-dump SRR14741309
fasterq-dump SRR14741310
fasterq-dump SRR14741311
fasterq-dump SRR14741312
fasterq-dump SRR14741313
fasterq-dump SRR14741314
fasterq-dump SRR14741315
fasterq-dump SRR14741316
fasterq-dump SRR14741317
fasterq-dump SRR14741318
fasterq-dump SRR14741319
fasterq-dump SRR14741320
fasterq-dump SRR14741321
fasterq-dump SRR14741322
fasterq-dump SRR14741323
fasterq-dump SRR14741324
fasterq-dump SRR14741325
fasterq-dump SRR14741326
fasterq-dump SRR14741327
fasterq-dump SRR14741328
fasterq-dump SRR14741329
fasterq-dump SRR14741330
fasterq-dump SRR14741331
fasterq-dump SRR14741332
fasterq-dump SRR14741333
fasterq-dump SRR14741334
fasterq-dump SRR14741335
fasterq-dump SRR14741336
fasterq-dump SRR14741337
fasterq-dump SRR14741338
fasterq-dump SRR14741339
fasterq-dump SRR14741340
fasterq-dump SRR14741341
fasterq-dump SRR14741342
fasterq-dump SRR14741343
fasterq-dump SRR14741344
fasterq-dump SRR14741345
fasterq-dump SRR14741346
fasterq-dump SRR14741347
fasterq-dump SRR14741348
fasterq-dump SRR14741349
fasterq-dump SRR14741350
fasterq-dump SRR14741351
fasterq-dump SRR14741352
fasterq-dump SRR14741353
fasterq-dump SRR14741354
fasterq-dump SRR14741355
fasterq-dump SRR14741356
fasterq-dump SRR14741357
fasterq-dump SRR14741358
fasterq-dump SRR14741359
fasterq-dump SRR14741360
fasterq-dump SRR14741361
fasterq-dump SRR14741362
fasterq-dump SRR14741363
fasterq-dump SRR14741364
fasterq-dump SRR14741365
fasterq-dump SRR14741366
fasterq-dump SRR14741367
fasterq-dump SRR14741368
fasterq-dump SRR14741369
fasterq-dump SRR14741370
fasterq-dump SRR14741371
fasterq-dump SRR14741372
fasterq-dump SRR14741373
fasterq-dump SRR14741374
fasterq-dump SRR14741375
fasterq-dump SRR14741376
fasterq-dump SRR14741377
fasterq-dump SRR14741378
fasterq-dump SRR14741379
fasterq-dump SRR14741380
fasterq-dump SRR14741381
fasterq-dump SRR14741382
fasterq-dump SRR14741383
fasterq-dump SRR14741384
fasterq-dump SRR14741385
fasterq-dump SRR14741386
fasterq-dump SRR14741387
fasterq-dump SRR14741388
fasterq-dump SRR14741389
fasterq-dump SRR14741390
fasterq-dump SRR14741391
fasterq-dump SRR14741392
fasterq-dump SRR14741393
fasterq-dump SRR14741394
fasterq-dump SRR14741395
fasterq-dump SRR14741396
fasterq-dump SRR14741397
fasterq-dump SRR14741398
fasterq-dump SRR14741399
fasterq-dump SRR14741400
fasterq-dump SRR14741401
fasterq-dump SRR14741402
fasterq-dump SRR14741403
fasterq-dump SRR14741404
fasterq-dump SRR14741405
fasterq-dump SRR14741406
fasterq-dump SRR14741407
fasterq-dump SRR14741408
fasterq-dump SRR14741409
fasterq-dump SRR14741410
fasterq-dump SRR14741411
fasterq-dump SRR14741412
fasterq-dump SRR14741413
fasterq-dump SRR14741414
fasterq-dump SRR14741415
fasterq-dump SRR14741416
fasterq-dump SRR14741417
fasterq-dump SRR14741418
fasterq-dump SRR14741419
fasterq-dump SRR14741420
fasterq-dump SRR14741421
fasterq-dump SRR14741422
fasterq-dump SRR14741423
fasterq-dump SRR14741424
fasterq-dump SRR14741425
fasterq-dump SRR14741426
fasterq-dump SRR14741427
fasterq-dump SRR14741428
fasterq-dump SRR14741429
fasterq-dump SRR14741430
fasterq-dump SRR14741431
fasterq-dump SRR14741432
fasterq-dump SRR14741433
fasterq-dump SRR14741434
fasterq-dump SRR14741435
fasterq-dump SRR14741436
fasterq-dump SRR14741437
fasterq-dump SRR14741438
fasterq-dump SRR14741439
fasterq-dump SRR14741440
fasterq-dump SRR14741441
fasterq-dump SRR14741442
fasterq-dump SRR14741443
fasterq-dump SRR14741444
fasterq-dump SRR14741445
fasterq-dump SRR14741446
fasterq-dump SRR14741447
fasterq-dump SRR14741448
fasterq-dump SRR14741449
fasterq-dump SRR14741450
fasterq-dump SRR14741451
fasterq-dump SRR14741452
fasterq-dump SRR14741453
fasterq-dump SRR14741454
fasterq-dump SRR14741455
fasterq-dump SRR14741456
fasterq-dump SRR14741457
fasterq-dump SRR14741458
fasterq-dump SRR14741459
fasterq-dump SRR14741460
fasterq-dump SRR14741461
fasterq-dump SRR14741462
fasterq-dump SRR14741463
fasterq-dump SRR14741464
fasterq-dump SRR14741465
fasterq-dump SRR14741466
fasterq-dump SRR14741467

echo "Please move the raw files into a suitable directory for running the R scripts."

# Run R scripts
Rscript ./R/01a_Process_Raw_Seqs_Run12.R
Rscript ./R/01b_Process_Raw_Seqs_Run13.R
Rscript ./R/01c_Process_Raw_Seqs_Run6.R
Rscript ./R/01d_Process_Raw_Seqs_Run7.R
Rscript ./R/02_Merge_and_Clean_Phyloseq_Objects.R
Rscript ./R/03_Add_Phylogeny.R
Rscript ./R/04_Explore_Phyloseq_Object.R
Rscript ./R/05_Beta_Diversity.R
Rscript ./R/06_Alpha_Diversity.R
Rscript ./R/07_Differential_Abundance.R
Rscript ./R/08_ml_models.R
Rscript ./R/08_Track_Read_Counts.R
Rscript ./R/09_Clean_Lee_2019-2020_metadata.R
Rscript ./R/10_Run_ITSxpress_on_Lee_2019-2020_Data.R
Rscript ./R/11_Process_Lee_2019-2020_ITS_raw_reads.R
Rscript ./R/12_combine_fungi_and_bacteria.R
Rscript ./R/13_compare_bact_and_fungi.R
Rscript ./R/14_Core_Microbiome.R
