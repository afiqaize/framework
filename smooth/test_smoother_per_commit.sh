#!/bin/bash

ftag="$1"
commonopt=('--tree' 'Events' '--file' '/nfs/dust/cms/group/exotica-desy/HeavyHiggs/Trees_mtt_B_new2/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8.root' '--variable' 'mtt : 320, 360, 400, 440, 480, 520, 560, 600, 640, 680, 720, 760, 800, 845, 890, 935, 985, 1050, 1140, 1300, 1700' '--variable' 'cHel = negate(b1k * b2k) + negate(b1r * b2r) + negate(b1n * b2n) : 5; -1, 1' '--weight' 'weight' '--bandwidth 0.2' '--bandwidth' '0.2')

if [ -z ${ftag} ]; then
    ftag='dev'
else
    git checkout "${ftag}"
fi

make -B smoother && ./smoother "${commonopt[@]}" --type weight --mode smooth --systematic TT_uF --weight 'wgtu = weight * MEfac_up' --weight 'wgtd = weight * MEfac_down' --output oweight_${ftag}.root && ./smoother "${commonopt[@]}" --type tree --mode smooth --systematic Jer --output otree_${ftag}.root && ./smoother "${commonopt[@]}" --mode histogram --output ohistogram_${ftag}.root
roothcat ohistogram_${ftag}.root mtt_cHel_unroll | grep sum > sum_histogram_${ftag}

for ifile in oweight_${ftag}.root otree_${ftag}.root; do
    roothcat ${ifile} mtt_cHel_unroll_nominal_source_template | grep sum >> sum_smooth_nominal_${ftag}
done

for ifile in oweight_${ftag}.root otree_${ftag}.root; do
    hname=$(rootfls ${ifile} | grep mtt_cHel_unroll_ | grep _up_smooth_template | awk '{print $3}' | awk -F';' '{print $1}')
    roothcat ${ifile} ${hname} | grep sum >> sum_smooth_up_${ftag}
done

echo ${ftag}
echo 'nominal histogram'
cat sum_histogram_${ftag}
rm sum_histogram_${ftag}
echo 'nominal smooth'
cat sum_smooth_nominal_${ftag}
rm sum_smooth_nominal_${ftag}
echo 'up smooth'
cat sum_smooth_up_${ftag}
rm sum_smooth_up_${ftag}
echo
echo
