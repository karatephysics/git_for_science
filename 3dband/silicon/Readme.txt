・main_brilian_igor.f90
2017/8/29作成
igorで3dplotをするために作成したもの。
DOSpara(1)=101のばあいは
igor上で
･Make/N=(100,100)/D w1
･w1='nan'
･SetScale/P y 0,0.01,"", w1
･SetScale/P x 0,0.01,"", w1
･w1=VB_SO[p+q*101] ただしVB_SOはロードしたエネルギーデータのファイル名(z軸のみ)
あとはimage plot、もしくはsurface plotで再現
詳しくは、/Users/Tetsu/Study/Optical_physics/Calculation/ARPES/silicon_3dband.pxp
参照
