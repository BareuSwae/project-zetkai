/* 프로젝트 젯카이 — 자동 채보: 오디오를 분석해 BPM·첫 박·난이도별 노트를 만든다 (개발자 도구 전용) */
(()=>{
'use strict';
const SR=22050,N=1024,HOP=256,FPS=SR/HOP,CTR=N*.75/SR*1000;

async function decode(blob){const Ctx=window.OfflineAudioContext||window.webkitOfflineAudioContext;const ac=new Ctx(1,1,SR);
 const buf=await ac.decodeAudioData(await blob.arrayBuffer());const n=buf.length,out=new Float32Array(n);
 for(let c=0;c<buf.numberOfChannels;c++){const d=buf.getChannelData(c);for(let i=0;i<n;i++)out[i]+=d[i]/buf.numberOfChannels}
 return{x:out,sr:buf.sampleRate}}

function fft(re,im){const n=re.length;for(let i=1,j=0;i<n;i++){let b=n>>1;for(;j&b;b>>=1)j^=b;j^=b;if(i<j){let t=re[i];re[i]=re[j];re[j]=t;t=im[i];im[i]=im[j];im[j]=t}}
 for(let len=2;len<=n;len<<=1){const a=-2*Math.PI/len,wr=Math.cos(a),wi=Math.sin(a);for(let i=0;i<n;i+=len){let cr=1,ci=0;for(let k=0;k<len/2;k++){const p=i+k,q=p+len/2,tr=re[q]*cr-im[q]*ci,ti=re[q]*ci+im[q]*cr;re[q]=re[p]-tr;im[q]=im[p]-ti;re[p]+=tr;im[p]+=ti;const nr=cr*wr-ci*wi;ci=cr*wi+ci*wr;cr=nr}}}}

/* onset envelopes (spectral flux) per band + RMS */
async function features(x,sr,prog){const hz=sr/N,b1=Math.max(2,Math.round(150/hz)),b2=Math.round(2000/hz),half=N/2;const frames=Math.max(1,Math.floor((x.length-N)/HOP));
 const win=new Float32Array(N);for(let i=0;i<N;i++)win[i]=.5-.5*Math.cos(2*Math.PI*i/(N-1));
 const lo=new Float32Array(frames),mi=new Float32Array(frames),hi=new Float32Array(frames),rms=new Float32Array(frames);let prev=new Float32Array(half),cur=new Float32Array(half);const re=new Float32Array(N),im=new Float32Array(N);
 for(let f=0;f<frames;f++){const o=f*HOP;let e=0;for(let i=0;i<N;i++){const v=x[o+i]||0;re[i]=v*win[i];im[i]=0;e+=v*v}rms[f]=Math.sqrt(e/N);fft(re,im);
  let a=0,b=0,c=0;for(let k=1;k<half;k++){const m=Math.log1p(100*Math.hypot(re[k],im[k]));cur[k]=m;const d=m-prev[k];if(d>0){if(k<=b1)a+=d;else if(k<=b2)b+=d;else c+=d}}
  lo[f]=a/b1;mi[f]=b/(b2-b1);hi[f]=c/(half-b2);const t=prev;prev=cur;cur=t;
  if(f%2000===0){prog&&prog(.1+.5*f/frames,'소리 분석 중…');await new Promise(r=>setTimeout(r,0))}}
 const norm=a=>{let m=1e-9;for(const v of a)m=Math.max(m,v);for(let i=0;i<a.length;i++)a[i]/=m};norm(lo);norm(mi);norm(hi);
 const env=new Float32Array(frames);for(let i=0;i<frames;i++)env[i]=lo[i]*1.1+mi[i]+hi[i]*.7;norm(env);
 return{lo,mi,hi,env,rms,frames}}

function tempo(env){const n=env.length;let mean=0;for(const v of env)mean+=v;mean/=n;const e=env.map(v=>v-mean);
 const lagOf=bpm=>60*FPS/bpm,lmin=Math.floor(lagOf(200)),lmax=Math.ceil(lagOf(60));const ac=new Float32Array(lmax+2);
 for(let l=lmin;l<=lmax+1;l++){let s=0;for(let i=0;i+l<n;i++)s+=e[i]*e[i+l];ac[l]=s/(n-l)}
 let best=lmin,bs=-1e9;for(let l=lmin;l<=lmax;l++){const bpm=60*FPS/l,w=Math.exp(-.5*Math.pow(Math.log2(bpm/125)/.6,2));const v=ac[l]*w+.5*(ac[2*l]||0)*w;if(v>bs){bs=v;best=l}}
 const y0=ac[best-1],y1=ac[best],y2=ac[best+1],den=y0-2*y1+y2,lag=best+(den?.5*(y0-y2)/den:0);let bpm=60*FPS/lag;
 // refine with a comb search around the estimate (fine BPM resolution)
 const comb=b=>{const P=60*FPS/b;let best=0;for(let ph=0;ph<P;ph+=.5){let v=0;for(let t=ph;t<n;t+=P){const i=t|0,f=t-i;v+=(env[i]||0)*(1-f)+(env[i+1]||0)*f}if(v>best)best=v}return best/(n/P)};
 let rb=bpm,rv=-1;for(let b=bpm-1.5;b<=bpm+1.5;b+=.05){const v=comb(b);if(v>rv){rv=v;rb=b}}bpm=rb;
 if(Math.abs(bpm-Math.round(bpm))<.3)bpm=Math.round(bpm);else bpm=Math.round(bpm*10)/10;return bpm}

function phase(env,lo,bpm){const P=60*FPS/bpm,n=env.length;let best=0,bs=-1;for(let ph=0;ph<P;ph+=.25){let s=0;for(let t=ph;t<n;t+=P){const i=t|0,f=t-i;s+=(env[i]||0)*(1-f)+(env[i+1]||0)*f}if(s>bs){bs=s;best=ph}}
 let first=0;for(let i=0;i<n;i++)if(env[i]>.25){first=i;break}
 let k=Math.floor((first-best)/P);if(best+k*P<first-P*.5)k++;let start=best+Math.max(0,k)*P;
 while(start-P>=first-P*.5)start-=P;return start/FPS*1000+CTR}

function peaks(F){const{env,lo,mi,hi}=F,n=env.length,W=Math.round(FPS*.2),out=[];let s=0;const pre=new Float32Array(n+1);for(let i=0;i<n;i++){s+=env[i];pre[i+1]=s}
 for(let i=2;i<n-2;i++){const v=env[i];if(v<.06)continue;let ok=true;for(let k=-3;k<=3;k++)if(k&&env[i+k]>v){ok=false;break}if(!ok)continue;
  const a=Math.max(0,i-W),b=Math.min(n,i+W),m=(pre[b]-pre[a])/(b-a);if(v<m*1.25+.04)continue;const band=lo[i]>=mi[i]&&lo[i]>=hi[i]*1.1?0:hi[i]>mi[i]*1.15?2:1;out.push({f:i,t:(i*HOP+N*.75)/SR*1000,s:v,band})}
 return out}

let seed=1;const rnd=()=>(seed=seed*16807%2147483647)/2147483647;
const DCFG={easy:{div:1,nps:1.25,min:1,chord:0,ln:3},normal:{div:2,nps:2.6,min:.5,chord:.92,ln:2},hard:{div:4,nps:4.6,min:.25,chord:.8,ln:1.5}};
function makeChart(pk,F,bpm,off,dk,durMs){const c=DCFG[dk],beat=60000/bpm,cand=new Map();
 for(const p of pk){const b=(p.t-off)/beat;if(b<0)continue;const q=Math.round(b*c.div)/c.div;if(Math.abs(b-q)*beat>Math.min(90,beat/c.div*.45))continue;const ex=cand.get(q);if(!ex||ex.s<p.s)cand.set(q,{b:q,s:p.s,band:p.band,f:p.f})}
 const allc=[...cand.values()].sort((a,b)=>a.b-b.b);let list=[...cand.values()].sort((a,b)=>b.s-a.s);const target=Math.max(8,Math.round(c.nps*durMs/1000));const picked=[];
 const near=(b)=>picked.some(x=>Math.abs(x.b-b)<c.min-1e-6);
 for(const p of list){if(picked.length>=target)break;if(near(p.b))continue;picked.push(p)}
 picked.sort((a,b)=>a.b-b.b);
 const sorted=[...picked].map(p=>p.s).sort((a,b)=>a-b),chordS=c.chord?sorted[Math.floor(sorted.length*c.chord)]||1:9;
 seed=Math.round(bpm*100+off)%2147483646+1+dk.length;const notes=[];let last=-1,lastB=-9,side=0;
 for(let i=0;i<picked.length;i++){const p=picked[i],gap=p.b-lastB;let pref=p.band===0?[0,3]:p.band===2?[1,2]:[0,1,2,3];
  pref=pref.filter(l=>!(gap<1&&l===last));if(gap<=.5)pref=pref.filter(l=>(l<2)!==(last<2)||last<0);if(!pref.length)pref=[0,1,2,3].filter(l=>l!==last);
  let l=pref[Math.floor(rnd()*pref.length)];if(pref.length>1&&rnd()<.6){const want=side?pref.filter(x=>x>=2):pref.filter(x=>x<2);if(want.length)l=want[Math.floor(rnd()*want.length)]}side=l<2?1:0;
  const n={b:p.b,l};const next=picked[i+1],ng=next?next.b-p.b:4;
  if(ng>=c.ln){const f0=p.f,f1=Math.min(F.rms.length-1,f0+Math.round((ng*beat*.7)/1000*FPS));let m=0,k=0;for(let f=f0+2;f<f1;f++){m+=F.rms[f];k++}m=k?m/k:0;
   if(m>F.rms[f0]*.55&&F.rms[f0]>.02){const e=p.b+Math.max(1,Math.floor((ng-.5)*c.div)/c.div);if(e>p.b)n.e=e}}
  notes.push(n);
  // keep the other hand busy during a long note: taps on other lanes, on the music's accents where there are any
  if(n.e!=null){const step=dk==='easy'?2:dk==='normal'?1:.5,avail=[0,1,2,3].filter(x=>x!==l),opp=avail.filter(x=>(x<2)!==(l<2));let lb=p.b,lx=-1;
   for(let b=p.b+step;b<n.e-step*.5+1e-6;b+=step){const nb=allc.filter(q=>q.b>lb+c.min-1e-6&&q.b<n.e-.24&&Math.abs(q.b-b)<step/2+1e-6);
    const bb=nb.length?nb.reduce((a,q)=>q.s>a.s?q:a).b:b;if(bb<=lb+1e-6||bb>=n.e-.24)continue;
    let pool=(rnd()<.7&&opp.length?opp:avail).filter(x=>x!==lx);if(!pool.length)pool=avail;const x=pool[Math.floor(rnd()*pool.length)];notes.push({b:bb,l:x});lx=x;lb=bb}}
  if(p.s>=chordS&&Math.abs(p.b-Math.round(p.b))<1e-6&&n.e==null){const l2=3-l!==l?3-l:(l+1)%4;notes.push({b:p.b,l:l2})}
  last=l;lastB=p.b}
 notes.sort((a,b)=>a.b-b.b||a.l-b.l);
 const secs=durMs/1000,nps=notes.length/Math.max(1,secs),lv=dk==='easy'?Math.min(5,Math.max(1,Math.round(1+nps*1.2))):dk==='normal'?Math.min(8,Math.max(3,Math.round(2+nps*1.2))):Math.min(13,Math.max(5,Math.round(3+nps*1.3)));
 return{level:lv,notes}}

/* main entry: returns {bpm, offset, charts:{easy,normal,hard}} */
window.AutoChart=async(blob,opt={})=>{const prog=opt.progress||(()=>{});prog(.02,'오디오 읽는 중…');
 const{x,sr}=await decode(blob);prog(.1,'소리 분석 중…');const F=await features(x,sr,prog);
 prog(.65,'박자 찾는 중…');await new Promise(r=>setTimeout(r,0));
 const bpm=opt.bpm||tempo(F.env);let offset=opt.offset!=null?opt.offset:phase(F.env,F.lo,bpm);
 prog(.8,'노트 고르는 중…');const pk=peaks(F);
 if(opt.offset==null&&pk.length>8){const beat=60000/bpm,ss=[...pk].sort((a,b)=>b.s-a.s).slice(0,Math.max(8,pk.length>>1)),res=ss.map(p=>{let r=((p.t-offset)%beat+beat)%beat;if(r>beat/2)r-=beat;return r}).filter(r=>Math.abs(r)<beat/6).sort((a,b)=>a-b);if(res.length>4)offset+=res[res.length>>1]}
const durMs=pk.length?pk[pk.length-1].t-pk[0].t:F.frames/FPS*1000;
 const charts={};for(const k of (opt.diffs||['easy','normal','hard']))charts[k]=makeChart(pk,F,bpm,offset,k,durMs);
 prog(1,'완성');return{bpm,offset:Math.round(offset),charts,duration:x.length/sr}};
})();
