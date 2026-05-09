# Fermentation Benchmark Result Analysis

분석 대상: 업로드된 `summary.csv`, `timeseries.csv`, `experiment_conditions.csv`, `environment_rankings.csv`.

- 총 조건 수: 108
- 환경별 조건 수: beer_wort=36, bread_dough=27, malt_mash_wort=18, rice_wine=27

## 1. 환경별 최상위 조건

### beer_wort
- Best condition: `BEER_ale_like_ale_warm_standard_wort_low_o2`
- Score=54.200, ethanol=49.178 g/L, CO2=51.639 g/L, retained CO2 4h=0.007 g/L
- Residual fermentable sugar=0.000 g/L, total carbohydrate consumed=106.000 g/L, max stress=0.134, failure mode=`completed`

### bread_dough
- Best condition: `BREAD_baker_like_standard_sugar_30C`
- Score=7.106, ethanol=11.626 g/L, CO2=12.530 g/L, retained CO2 4h=2.689 g/L
- Residual fermentable sugar=3.921 g/L, total carbohydrate consumed=26.303 g/L, max stress=0.025, failure mode=`proof_sugar_depleted`

### malt_mash_wort
- Best condition: `MALT_stress_tolerant_standard_mash_62C_robust25_heavy_grist`
- Score=62.529, ethanol=38.050 g/L, CO2=39.885 g/L, retained CO2 4h=0.005 g/L
- Residual fermentable sugar=59.410 g/L, total carbohydrate consumed=97.961 g/L, max stress=0.415, failure mode=`nitrogen_limited`

### rice_wine
- Best condition: `RICE_stress_tolerant_strong_koji_25C`
- Score=83.059, ethanol=87.816 g/L, CO2=92.167 g/L, retained CO2 4h=0.022 g/L
- Residual fermentable sugar=6.685 g/L, total carbohydrate consumed=206.175 g/L, max stress=0.227, failure mode=`completed`


## 2. 핵심 해석

- 제빵환경은 `baker_like + 30°C`가 가장 강하다. 30°C 평균 objective가 가장 높고, 22°C는 CO₂ 생산이 뚜렷하게 느리다.
- 쌀술은 `strong_koji`와 `stress_tolerant` 조합이 가장 강하다. 25°C와 30°C는 비슷하지만 25°C top condition이 completion 판정을 받아 논문용 대표 조건으로 더 깔끔하다.
- 맥주는 `standard_wort + low_o2 + ale_like`가 가장 안정적이다. high gravity 조건은 대부분 nitrogen-limited로 분류되어, 추가 실험축으로 FAN/질소 보강을 넣기 좋다.
- 맥아/매싱 계열은 62°C mash가 압도적으로 좋고, 55°C/68°C는 starch conversion limited가 많다. 다만 최고 조건도 nitrogen-limited라서 다음 sweep에서는 protease/nitrogen 조건을 키우는 것이 좋다.


## 3. 지표 주의사항

- 쌀술과 맥아는 동적 당화가 포함되므로 `apparent_attenuation_fraction`, `fermentable_sugar_used_g_L`가 초기 발효당 기준에서 왜곡될 수 있다. 그래서 `corrected_summary.csv`에 total carbohydrate 기준 보정 지표를 추가했다.
- 논문 본문에서는 맥주처럼 이미 당화된 wort는 attenuation을 사용하고, 쌀술/맥아처럼 당화-발효 결합계는 `total_carbohydrate_conversion_fraction`, `ethanol_yield_per_total_carbohydrate_consumed_g_g`, `residual_total_carbohydrate_fraction`을 쓰는 편이 안전하다.


## 4. 실패 모드 분포

- beer_wort: completed=18, nitrogen_limited=12, balanced_or_slow=6
- bread_dough: proof_residual_sugar_expected=16, proof_sugar_depleted=11
- malt_mash_wort: mash_starch_conversion_limited=12, nitrogen_limited=5, balanced_or_slow=1
- rice_wine: saccharification_limited=14, completed=5, nitrogen_limited=4, balanced_or_slow=2, ethanol_inhibited=2

## 5. 다음 실험 제안

1. Bread: baker_like, 30°C 주변에서 26/28/30/32/34°C 미세 sweep + CO₂ retention coefficient sweep.
2. Rice wine: strong_koji, stress_tolerant, 24–30°C sweep + 초기 pH 3.8–4.8 + 질소원 보강.
3. Beer: standard/high-gravity wort에서 nitrogen_g_L/FAN proxy를 늘려 nitrogen-limited 해소 여부 확인.
4. Malt: 62°C 중심으로 mash time 1/2/4h, protease activity, initial nitrogen을 추가 sweep.
