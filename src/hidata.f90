Block Data hidata
  Parameter (maxstr=150001)
  Double Precision xl(10), xu(10), acc
  Common /bveg1/xl, xu, acc, ndim, ncall, itmx, nprn
  Common /sedvax/num1
  Common /hparnt/hipr1(100), ihpr2(50), hint1(100), ihnt2(50)
  Common /hmain1/eatt, jatt, natt, nt, np, n0, n01, n10, n11
  Common /hmain2/katt(maxstr, 4), patt(maxstr, 4)
  Common /hstrng/nfp(300, 15), pp(300, 15), nft(300, 15), pt(300, 15)
  Common /hjcrdn/yp(3, 300), yt(3, 300)
  Common /hjjet1/npj(300), kfpj(300, 500), pjpx(300, 500), pjpy(300, 500), pjpz(300, 500), pjpe(300, 500), pjpm(300, 500), ntj(300), kftj(300, 500), pjtx(300, 500), pjty(300, 500), pjtz(300, 500), pjte(300, 500), pjtm(300, 500)
  Common /hjjet2/nsg, njsg(maxstr), iasg(maxstr, 3), k1sg(maxstr, 100), k2sg(maxstr, 100), pxsg(maxstr, 100), pysg(maxstr, 100), pzsg(maxstr, 100), pesg(maxstr, 100), pmsg(maxstr, 100)
  Common /hijdat/hidat0(10, 10), hidat(10)
  Common /hpint/mint4, mint5, atco(200, 20), atxs(0:200)
  Save
  Data num1/30123984/, xl/10*0.D0/, xu/10*1.D0/
  Data ncall/1000/, itmx/100/, acc/0.01/, nprn/0/
  Data hipr1/1.5, 0.35, 0.5, 0.9, 2.0, 0.1, 1.5, 2.0, -1.0, -2.25, 2.0, 0.5, 1.0, 2.0, 0.2, 2.0, 2.5, 0.3, 0.1, 1.4, 1.6, 1.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.4, 57.0, 28.5, 3.9, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.14159, 0.0, 0.4, 0.1, 1.5, 0.1, 0.25, 0.0, 0.5, 0.0, 0.0, 50*0.0/
  Data ihpr2/1, 3, 0, 1, 1, 1, 1, 10, 0, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 30*0/
  Data hint1/100*0/
  Data ihnt2/50*0/
  Data natt/0/, eatt/0.0/, jatt/0/, nt/0/, np/0/, n0/0/, n01/0/, n10/0/, n11/0/
  Data katt/600004*0/, patt/600004*0.0/
  Data nfp/4500*0/, pp/4500*0.0/, nft/4500*0/, pt/4500*0.0/
  Data yp/900*0.0/, yt/900*0.0/
  Data npj/300*0/, kfpj/150000*0/, pjpx/150000*0.0/, pjpy/150000*0.0/, pjpz/150000*0.0/, pjpe/150000*0.0/, pjpm/150000*0.0/
  Data ntj/300*0/, kftj/150000*0/, pjtx/150000*0.0/, pjty/150000*0.0/, pjtz/150000*0.0/, pjte/150000*0.0/, pjtm/150000*0.0/
  Data nsg/0/, njsg/150001*0/, iasg/450003*0/, k1sg/15000100*0/, k2sg/15000100*0/, pxsg/15000100*0.0/, pysg/15000100*0.0/, pzsg/15000100*0.0/, pesg/15000100*0.0/, pmsg/15000100*0.0/
  Data mint4/0/, mint5/0/, atco/4000*0.0/, atxs/201*0.0/
  Data (hidat0(1,i), i=1, 10)/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.25, 2.5, 4.0, 4.1/
  Data (hidat0(2,i), i=1, 10)/2.0, 3.0, 5.0, 6.0, 7.0, 8.0, 8.0, 10.0, 10.0, 10.0/
  Data (hidat0(3,i), i=1, 10)/1.0, 0.8, 0.8, 0.7, 0.45, 0.215, 0.21, 0.19, 0.19, 0.19/
  Data (hidat0(4,i), i=1, 10)/0.35, 0.35, 0.3, 0.3, 0.3, 0.3, 0.5, 0.6, 0.6, 0.6/
  Data (hidat0(5,i), i=1, 10)/23.8, 24.0, 26.0, 26.2, 27.0, 28.5, 28.5, 28.5, 28.5, 28.5/
  Data ((hidat0(j,i),i=1,10), j=6, 9)/40*0.0/
  Data (hidat0(10,i), i=1, 10)/5.0, 20.0, 53.0, 62.0, 100.0, 200.0, 546.0, 900.0, 1800.0, 4000.0/
  Data hidat/10*0.0/
End Block Data hidata
