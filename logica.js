// Utilidades de números y complejos
    const EPS = 1e-10;
    const TOL_CLUSTER = 1e-6; // para agrupar raíces por cercanía

    // nuevo: formatea escalares: si es casi entero -> entero; else recorta decimales; usa exponencial si aplica
    function fmtScalar(x, dp=6){
      const ax = Math.abs(x);
      if (!Number.isFinite(x)) return String(x);
      // entero cercano
      if (Math.abs(x - Math.round(x)) < 1e-12) return String(Math.round(x));
      if (ax < 1e-12) return '0';
      if (ax >= 1e6 || ax < 1e-4) return x.toExponential(6);
      // recorta a dp decimales y elimina ceros/trailing dot
      let s = x.toFixed(dp);
      // quitar ceros finales y posible punto
      s = s.replace(/\.?0+$/,'')
      return s;
    }

    // nuevo: devuelve la expresión que representa "coef * x" sin mostrar 1: 1 -> "x", -1 -> "-x", otros -> "k x"
    function coefTimesX(a){
      if (!Number.isFinite(a)) return fmtNum(a) + ' x';
      if (Math.abs(a) < 1e-12) return '0';
      if (Math.abs(a - 1) < 1e-12) return 'x';
      if (Math.abs(a + 1) < 1e-12) return '-x';
      return `${fmtNum(a)} x`;
    }

    // Reemplazar coefTimesX para salida HTML (ej: "3x", "x", "-x")
    function coefTimesXHtml(a){
      if (!Number.isFinite(a)) return fmtNum(a) + 'x';
      if (Math.abs(a) < 1e-12) return '0';
      if (Math.abs(a - 1) < 1e-12) return 'x';
      if (Math.abs(a + 1) < 1e-12) return '-x';
      return `${fmtNum(a)}x`;
    }

    class C {
      constructor(re=0, im=0){ this.re = re; this.im = im; }
      static from(x){ return x instanceof C ? x : new C(x,0); }
      add(z){ return new C(this.re+z.re, this.im+z.im); }
      sub(z){ return new C(this.re-z.re, this.im-z.im); }
      mul(z){ return new C(this.re*z.re - this.im*z.im, this.re*z.im + this.im*z.re); }
      div(z){
        const d = z.re*z.re + z.im*z.im; return new C((this.re*z.re + this.im*z.im)/d, (this.im*z.re - this.re*z.im)/d);
      }
      abs(){ return Math.hypot(this.re, this.im); }
      eq(z, tol=TOL_CLUSTER){ return Math.hypot(this.re - z.re, this.im - z.im) < tol; }
      toString(dp=6){
        // usa fmtScalar para mostrar enteros sin decimales
        const r = Math.abs(this.re) < 1e-12 ? 0 : this.re;
        const i = Math.abs(this.im) < 1e-12 ? 0 : this.im;
        if (Math.abs(i) < 1e-12) return fmtScalar(r, dp);
        if (Math.abs(r) < 1e-12) return (i>=0?"":"-") + 'i' + fmtScalar(Math.abs(i), dp);
        const sign = i>=0? "+": "-";
        return `${fmtScalar(r, dp)} ${sign} i${fmtScalar(Math.abs(i), dp)}`;
      }
    }

    // Evalúa polinomio p(r) = a0*r^n + a1*r^(n-1) + ... + an
    function polyEval(coeffs, r){
      // coeffs: [a_n, a_{n-1}, ..., a_0]
      r = C.from(r);
      let x = new C(0,0);
      for (let i=0; i<coeffs.length; i++){
        x = x.mul(r).add(new C(coeffs[i],0));
      }
      return x;
    }

    // Durand–Kerner para raíces de polinomio (monic). Devuelve array de C
    function durandKerner(coeffs, maxIter=200, tol=1e-12){
      // Normaliza a monic: p(z) = z^n + b_{n-1} z^{n-1} + ... + b0
      const a = coeffs.slice();
      if (Math.abs(a[0]) < EPS) throw new Error("El coeficiente líder no puede ser 0.");
      for (let i=1;i<a.length;i++) a[i] /= a[0];
      a[0] = 1;
      const n = a.length-1;

      // p(z) como monic
      const p = (z)=>{
        z = C.from(z);
        let y = new C(1,0); // z^0
        let s = new C(1,0); // comenzaremos desde el término más bajo si queremos, pero usaremos Horner inverso
        // mejor usar Horner de mayor a menor: a0*z^n + a1*z^{n-1} ...
        let acc = new C(0,0);
        for (let i=0;i<=n;i++){
          acc = acc.mul(z).add(new C(a[i],0));
        }
        return acc; // monic
      };

      // Inicialización: raíces en círculo
      const roots = [];
      const R = 1; // radio inicial
      for (let k=0;k<n;k++){
        const angle = 2*Math.PI*k/n;
        roots.push(new C(R*Math.cos(angle), R*Math.sin(angle)));
      }

      for (let it=0; it<maxIter; it++){
        let maxDelta = 0;
        for (let i=0;i<n;i++){
          let denom = new C(1,0);
          for (let j=0;j<n;j++) if (j!==i) denom = denom.mul(roots[i].sub(roots[j]));
          const pi = p(roots[i]);
          const delta = pi.div(denom);
          roots[i] = roots[i].sub(delta);
          maxDelta = Math.max(maxDelta, delta.abs());
        }
        if (maxDelta < tol) break;
      }
      return roots;
    }

    function parseCoeffs(){
      const raw = document.getElementById('coeffs').value.trim();
      if (!raw) throw new Error('Ingresa los coeficientes separados por comas.');
      const parts = raw.split(/[,\s]+/).map(Number);
      if (parts.some(x=>!Number.isFinite(x))) throw new Error('Hay elementos no numéricos.');
      if (parts.length < 2) throw new Error('Se requieren al menos dos coeficientes (orden ≥ 1).');
      if (Math.abs(parts[0]) < EPS) throw new Error('El primer coeficiente (a_n) no puede ser 0.');
      return parts;
    }

    // Reemplazar charPolyString para devolver HTML con potencias en <sup>
    function charPolyStringHtml(coeffs){
      const n = coeffs.length - 1;
      const terms = [];
      for (let i=0;i<coeffs.length;i++){
        const a = coeffs[i];
        const pow = n - i;
        if (Math.abs(a) < EPS) continue;
        let aStr;
        if (Math.abs(a) === 1 && pow > 0) aStr = (a>0 ? "" : "-");
        else aStr = fmtScalar(a);
        const rStr = pow===0? "" : (pow===1? "r" : `r<sup>${pow}</sup>`);
        const term = `${aStr}${rStr}`.trim();
        terms.push(term);
      }
      const s = (terms.length? terms.join(' + ').replace(/\+ -/g,' - ') : '0') + ' = 0';
      return `<div class="text-slate-800"><code class="font-mono">${s}</code></div>`;
    }

    // Actualiza generalSolution para producir HTML bonito con sub/sup y C<sub>i</sub>
    function generalSolutionHtml(clusters){
      const reals = clusters.filter(c=>c.type==='real').sort((a,b)=>b.mult-a.mult || a.root.re-b.root.re);
      const complexes = clusters.filter(c=>c.type==='complex');

      const used = new Array(complexes.length).fill(false);
      const pairs = [];
      for (let i=0;i<complexes.length;i++){
        if (used[i]) continue;
        const a = complexes[i];
        let jfound = -1;
        for (let j=i+1;j<complexes.length;j++){
          if (used[j]) continue;
          const b = complexes[j];
          if (Math.abs(a.root.re - b.root.re) < TOL_CLUSTER && Math.abs(a.root.im + b.root.im) < TOL_CLUSTER){
            jfound = j; break;
          }
        }
        if (jfound === -1){
          pairs.push({ alpha: a.root.re, beta: Math.abs(a.root.im), mult: a.mult });
          used[i] = true;
        } else {
          const b = complexes[jfound];
          pairs.push({ alpha: (a.root.re + b.root.re)/2, beta: (Math.abs(a.root.im) + Math.abs(b.root.im))/2, mult: Math.min(a.mult, b.mult) });
          used[i] = used[jfound] = true;
        }
      }

      const parts = [];
      let cidx = 1;
      // Reales
      for (const c of reals){
        for (let k=0;k<c.mult;k++){
          const Cn = `C<sub>${cidx++}</sub>`;
          const xpow = k===0? '' : (k===1? 'x' : `x<sup>${k}</sup>`);
          const alpha = c.root.re;
          if (Math.abs(alpha) < 1e-12){
            parts.push(`<span class="inline-block">${Cn} ${xpow}</span>`);
          } else {
            // usar fracTimesXHtml para evitar decimales largos; devuelve p.ej. "3/2x"
            const expArg = fracTimesXHtml(alpha);
            parts.push(`<span class="inline-block">${Cn} ${xpow}e<sup>${expArg}</sup></span>`);
          }
        }
      }

      // Complejas por pares
      for (const p of pairs){
        for (let k=0;k<p.mult;k++){
          const Cn = `C<sub>${cidx++}</sub>`;
          const Dn = `C<sub>${cidx++}</sub>`;
          const xpow = k===0? '' : (k===1? 'x' : `x<sup>${k}</sup>`);
          // usar fracTimesXHtml para beta y alpha
          const betaExpr = fracTimesXHtml(p.beta);
          if (Math.abs(p.alpha) < 1e-12){
            parts.push(`<span class="inline-block">${xpow}(${Cn} cos(${betaExpr}) + ${Dn} sin(${betaExpr}))</span>`);
          } else {
            const alphaExpr = fracTimesXHtml(p.alpha);
            parts.push(`<span class="inline-block">e<sup>${alphaExpr}</sup>${xpow}(${Cn} cos(${betaExpr}) + ${Dn} sin(${betaExpr}))</span>`);
          }
        }
      }

      if (parts.length===0) return `<div class="text-slate-700">—</div>`;
      // unir con + y mostrar en bloque <pre> para preservar espaciamiento monoespaciado
      return `<div class="text-slate-800"><div class="mb-2 font-medium">y(x) =</div><div class="flex flex-wrap gap-3">${parts.map(p=>`<div class="chip">${p}</div>`).join(' + ')}</div></div>`;
    }

    // Agrupa raíces por cercanía (para multiplicidades). Devuelve [{root: C, mult: k, type: 'real'|'complex'}]
    function clusterRoots(roots){
      const clusters = [];
      roots.forEach(r => {
        // intenta unir a un cluster existente
        let found = false;
        for (const c of clusters){
          if (r.abs()===Infinity) continue;
          if (r.eq(c.root)){
            c.items.push(r);
            // actualiza promedio
            const m = c.items.length;
            const sum = c.items.reduce((acc,z)=> new C(acc.re+z.re, acc.im+z.im), new C(0,0));
            c.root = new C(sum.re/m, sum.im/m);
            found = true; break;
          }
        }
        if (!found) clusters.push({root: new C(r.re, r.im), items: [r]});
      });
      // mapear
      return clusters.map(c=>({ root: c.root, mult: c.items.length, type: Math.abs(c.root.im) < 1e-8 ? 'real' : 'complex' }));
    }

    // Genera string de solución general
    function generalSolution(clusters){
      // Orden: reales (de mayor multiplicidad a menor), luego complejas (por par)
      const reals = clusters.filter(c=>c.type==='real').sort((a,b)=>b.mult-a.mult || a.root.re-b.root.re);
      const complexes = clusters.filter(c=>c.type==='complex');

      // Para complejas, emparejar conjugados y sumar multiplicidades (asumiendo pares correctos)
      const used = new Array(complexes.length).fill(false);
      const pairs = [];
      for (let i=0;i<complexes.length;i++){
        if (used[i]) continue;
        const a = complexes[i];
        // busca conjugado
        let jfound = -1;
        for (let j=i+1;j<complexes.length;j++){
          if (used[j]) continue;
          const b = complexes[j];
          if (Math.abs(a.root.re - b.root.re) < TOL_CLUSTER && Math.abs(a.root.im + b.root.im) < TOL_CLUSTER){
            jfound = j; break;
          }
        }
        if (jfound === -1){
          // si no se encuentra, lo tratamos solo (puede ser error numérico); lo incluimos como par con su conjugado implícito
          pairs.push({ alpha: a.root.re, beta: Math.abs(a.root.im), mult: a.mult });
          used[i] = true;
        } else {
          const b = complexes[jfound];
          pairs.push({ alpha: (a.root.re + b.root.re)/2, beta: (Math.abs(a.root.im) + Math.abs(b.root.im))/2, mult: Math.min(a.mult, b.mult) });
          used[i] = used[jfound] = true;
        }
      }

      let parts = [];
      let cidx = 1;

      // Reales
      for (const c of reals){
        for (let k=0;k<c.mult;k++){
          const Cn = `C_${cidx++}`;
          const xpow = k===0? '' : (k===1? 'x ' : `x^${k} `);
          const alpha = c.root.re;
          if (Math.abs(alpha) < 1e-12){
            // si exponente ~0, omitimos e^{0 x}
            parts.push(`${(Cn)} ${xpow}`.trim());
          } else {
            // usar fracTimesX en lugar de coefTimesX para fracciones
            const expArg = fracTimesX(alpha); // por ejemplo "3/2 x"
            parts.push(`${Cn} ${xpow}e^{${expArg}}`.trim());
          }
        }
      }

      // Complejas por pares (con multiplicidad)
      for (const p of pairs){
        for (let k=0;k<p.mult;k++){
          const Cn = `C_${cidx++}`;
          const Dn = `C_${cidx++}`;
          const xpow = k===0? '' : (k===1? 'x ' : `x^${k} `);
          const betaExpr = fracTimesX(p.beta); // "k x" o fracción
          if (Math.abs(p.alpha) < 1e-12){
            parts.push(`${xpow}(${Cn} \\cos(${betaExpr}) + ${Dn} \\sin(${betaExpr}))`.trim());
          } else {
            const alphaExpr = fracTimesX(p.alpha);
            parts.push(`e^{${alphaExpr}}${xpow}(${Cn} \\cos(${betaExpr}) + ${Dn} \\sin(${betaExpr}))`.trim());
          }
        }
      }

      if (parts.length===0) return '—';
      return 'y(x) = ' + parts.join(' + ');
    }

// Añadir utilitario fmtNum (usa fmtScalar existente)
function fmtNum(x){ return fmtScalar(x, 6); }

// Selecciones: escribir dentro del contenedor .out para no romper la estructura del acordeón
const elChar = document.getElementById('charpoly'); // ya es el .out interno
const elRoots = document.querySelector('#roots .out'); // Cambiado: apuntar al .out interno
const elSol = document.querySelector('#solution .out'); // Cambiado: apuntar al .out interno

// --- Nuevo: sistema de notificaciones (toasts) ---
function ensureToastContainer() {
  let c = document.getElementById('toast-container');
  if (c) return c;
  c = document.createElement('div');
  c.id = 'toast-container';
  Object.assign(c.style, {
    position: 'fixed',
    right: '20px',
    bottom: '24px',
    display: 'flex',
    flexDirection: 'column',
    gap: '10px',
    alignItems: 'flex-end',
    zIndex: 9999,
    pointerEvents: 'none'
  });
  document.body.appendChild(c);
  return c;
}

function showToast(message, type='error', duration=4500) {
  const container = ensureToastContainer();
  const t = document.createElement('div');
  t.className = 'toast';
  // estilos inline para evitar depender de CSS externo
  Object.assign(t.style, {
    pointerEvents: 'auto',
    minWidth: '220px',
    maxWidth: '420px',
    padding: '10px 14px',
    borderRadius: '10px',
    color: '#061220',
    fontWeight: 600,
    boxShadow: '0 8px 24px rgba(2,8,23,0.6)',
    transform: 'translateY(8px)',
    opacity: '0',
    transition: 'opacity 260ms ease, transform 260ms cubic-bezier(.2,.9,.2,1)'
  });

  // color según tipo
  if (type === 'success') {
    t.style.background = 'linear-gradient(90deg, #10b981, #059669)';
    t.style.color = '#e6fffa';
  } else if (type === 'warn') {
    t.style.background = 'linear-gradient(90deg, #f59e0b, #f97316)';
    t.style.color = '#071422';
  } else { // error / default
    t.style.background = 'linear-gradient(90deg, #fb7185, #fb923c)';
    t.style.color = '#061220';
  }

  t.innerText = message;
  container.appendChild(t);

  // aparecer
  requestAnimationFrame(() => {
    t.style.opacity = '1';
    t.style.transform = 'translateY(0)';
  });

  // auto dismiss
  const hide = () => {
    t.style.opacity = '0';
    t.style.transform = 'translateY(8px)';
    setTimeout(()=> t.remove(), 300);
  };
  const timeoutId = setTimeout(hide, duration);

  // cerrar al hacer click
  t.addEventListener('click', () => {
    clearTimeout(timeoutId);
    hide();
  });
}

// Reemplazar showErrorHtml para usar toasts en lugar de escribir en el panel
function showErrorHtml(el, msg){
  // Mostrar notificación flotante
  showToast(msg, 'warn', 5000);
  // Opcional: limpiar/ocultar contenido del panel afectado
  try {
    if (el && el instanceof HTMLElement) el.innerHTML = '—';
  } catch(e) { /* silencioso */ }
}

// Reemplazar mostrarChar/calcularRaices/mostrarSolucion para inyectar HTML
    function mostrarChar(){
      try{
        const coeffs = parseCoeffs();
        elChar.innerHTML = charPolyStringHtml(coeffs);
      }catch(e){
        showErrorHtml(elChar, e.message);
      }
    }

    function calcularRaices(){
      try{
        const coeffs = parseCoeffs();
        const roots = durandKerner(coeffs);
        const items = roots.map((z,i)=>{
          const idx = i+1;
          const tag = Math.abs(z.im) < 1e-8 ? 'real' : 'compleja';
          // formatear usando fracciones cuando sea posible
          const val = formatComplexAsFraction(z);
          return `<li class="py-2 flex flex-col sm:flex-row sm:items-center sm:justify-between border-b last:border-b-0">
                    <div class="text-slate-800"><span class="font-semibold">r<sub>${idx}</sub></span> = <span class="chip font-mono">${val}</span></div>
                    <div class="text-sm text-slate-500 mt-1 sm:mt-0">${tag}</div>
                  </li>`;
        }).join('');
        elRoots.innerHTML = `<ol class="list-decimal pl-5 space-y-1 bg-slate-50/50 rounded p-2">${items}</ol>`;
      }catch(e){
        showErrorHtml(elRoots, e.message);
      }
    }

    function mostrarSolucion(){
      try{
        const coeffs = parseCoeffs();
        const roots = durandKerner(coeffs);
        const clusters = clusterRoots(roots);
        elSol.innerHTML = generalSolutionHtml(clusters);
      }catch(e){
        showErrorHtml(elSol, e.message);
      }
    }

    // Hacer funciones globales (si no están ya)
    window.mostrarChar = mostrarChar;
    window.calcularRaices = calcularRaices;
    window.mostrarSolucion = mostrarSolucion;

    // nuevo: aproximar un número real por una fracción usando fracciones continuas
function toFraction(x, maxDen=1000, tol=1e-10){
  if (!Number.isFinite(x)) return null;
  const sign = x < 0 ? -1 : 1;
  x = Math.abs(x);
  const a = [];
  let q = x;
  for (let i=0;i<20;i++){
    const ai = Math.floor(q);
    a.push(ai);
    const frac = a.reduceRight((acc, v, idx)=>{
      if (idx===a.length-1) return v;
      return v + 1/acc;
    }, null);
    // compute numerator/denominator from continued fraction
    // but easier: compute convergents
    // compute convergent p/q
    let p0=1, q0=0, p1=a[0], q1=1;
    for (let k=1;k<a.length;k++){
      const ak = a[k];
      const p2 = ak*p1 + p0;
      const q2 = ak*q1 + q0;
      p0 = p1; q0 = q1; p1 = p2; q1 = q2;
    }
    if (q1 > maxDen) break;
    const approx = p1/q1;
    if (Math.abs(approx - x) <= Math.max(tol, Math.abs(x)*1e-12)) {
      return {num: sign*p1, den: q1};
    }
    const rem = q - ai;
    if (rem === 0) break;
    q = 1/rem;
  }
  // fallback: try simple rationalization by searching denominators up to maxDen
  for (let den=1; den<=maxDen; den++){
    const num = Math.round(x*den);
    if (Math.abs(num/den - x) <= Math.max(tol, Math.abs(x)*1e-12)){
      return {num: sign*num, den};
    }
  }
  return null;
}

// nuevo: formatea número como fracción si posible, else decimal con fmtNum
function fracOrDecimal(x, maxDen=1000, tol=1e-8){
  if (!Number.isFinite(x)) return String(x);
  // considerar entero cercano
  if (Math.abs(x - Math.round(x)) < 1e-12) return String(Math.round(x));
  // intentar fracción directa
  const f = toFraction(x, maxDen, tol);
  if (f){
    if (f.den === 1) return String(f.num);
    return `${f.num}/${f.den}`;
  }
  // intentar representar x como sqrt(p/q) si x^2 es racional aproximable
  const sf = toFraction(x*x, maxDen, tol);
  if (sf){
    const fracStr = sf.den === 1 ? `${sf.num}` : `${sf.num}/${sf.den}`;
    // si x es negativo, mostrar -√(...)
    return (x < 0 ? '-' : '') + `√(${fracStr})`;
  }
  return fmtNum(x);
}

// nuevo: formatea número *x* (ej: "x", "-x", "3/2 x", "3x")
function fracTimesXHtml(a){
  if (!Number.isFinite(a)) return fmtNum(a) + 'x';
  if (Math.abs(a) < 1e-12) return '0';
  if (Math.abs(a - 1) < 1e-12) return 'x';
  if (Math.abs(a + 1) < 1e-12) return '-x';
  const s = fracOrDecimal(a, 200);
  return `${s}x`;
}

// nuevo: formatea número *a* para "k x" en texto (usa fracciones si aplica)
function fracTimesX(a){
  if (!Number.isFinite(a)) return fmtNum(a) + ' x';
  if (Math.abs(a) < 1e-12) return '0';
  if (Math.abs(a - 1) < 1e-12) return 'x';
  if (Math.abs(a + 1) < 1e-12) return '-x';
  const s = fracOrDecimal(a, 200); // usa aproximación a fracción
  return `${s} x`;
}

// nuevo: formatea complejo z (C) usando fracciones cuando aplica
function formatComplexAsFraction(z){
  const re = Math.abs(z.re) < 1e-12 ? 0 : z.re;
  const im = Math.abs(z.im) < 1e-12 ? 0 : z.im;
  if (im === 0) return fracOrDecimal(re, 200);
  const imAbsStr = fracOrDecimal(Math.abs(im), 200);
  // Formato para parte imaginaria: "3/2i"
  const imagPart = `${imAbsStr}i`;
  if (re === 0) {
    return (im > 0) ? imagPart : `-${imagPart}`;
  }
  const reStr = fracOrDecimal(re, 200);
  const sign = im >= 0 ? ' + ' : ' - ';
  return `${reStr}${sign}${imagPart}`;
}