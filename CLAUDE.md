# CranioPlan — reglas del proyecto para Claude Code

Tesis de Ingeniería Biomédica (FCEFyN, UNC). Autores: Valentino Andri e Ignacio
Gross (Nacho). Colaborador clínico: Servicio de Neurocirugía, Hospital Garrahan.
Módulo de 3D Slicer 5.10 para planificación prequirúrgica de craneosinostosis.

## Arquitectura del módulo (paquete, no archivo único)

- `CranioPlan.py`: módulo de Slicer con tres clases. `CranioPlan` (registro),
  `CranioPlanWidget` (interfaz, sin procesamiento) y `CranioPlanLogic`
  (orquestación, sin elementos de interfaz).
- `CranioPlanLib/`: un archivo por bloque (`BloqueC.py` carga DICOM,
  `BloqueA.py` cráneo, `BloqueF.py` corte, `BloqueG.py` reacomodamiento),
  `Comun.py` con nombres de nodos y logger, `__init__.py` con `recargar()`.
- Los bloques se comunican por NOMBRE DE NODO en la escena, definidos solo en
  `Comun.py`. No inventar nombres nuevos fuera de ahí.
- Cada bloque debe seguir siendo pegable tal cual en la consola de Python de
  Slicer (por eso las constantes quedan como variables de módulo).
- Si se agrega un archivo a `CranioPlanLib/`, agregarlo también en
  `CMakeLists.txt`.
- El núcleo de decisión del Bloque A es el mismo código que el script de
  diagnóstico del relevamiento: no divergirlos.

## Estilo de trabajo

- Código, comentarios y textos de interfaz en español, sin jerga para el médico.
- Explicar la decisión de arquitectura ANTES de implementar cambios grandes.
- Documentar en el código los casos conocidos en que un algoritmo falla.
- Preferir pipelines vectorizados VTK/NumPy; evitar bucles Python sobre
  vértices o vóxeles.
- Slicer 5.10: los slots de Qt no aceptan kwargs (`addWidget(w, 2)`, no
  `stretch=2`).
- Versionado vX.Y por bloque; registrar cada cambio en `CHANGELOG.txt`.
- Si algo no tiene sentido técnico o hay una opción mejor, decirlo antes de
  implementar.

## Git (repo compartido con Nacho)

- Nunca commitear directo a `master` ni a `main`: rama por tarea + Pull Request.
- Nunca `git push --force`. Nunca reescribir historia ya subida.
- `git pull` al empezar cada sesión.
- Para mover o renombrar archivos usar `git mv` (conserva el historial).
- Commits chicos y con mensaje en español que diga qué y por qué.

## Datos de pacientes

PROHIBIDO leer, copiar, mover o commitear estudios DICOM, NRRD, NIfTI, MRB, ni
CSV/Excel con datos de casos. Son datos reales de pacientes pediátricos.

## Pruebas

Claude no puede ejecutar la interfaz de Slicer ni tiene los estudios. Las
pruebas end-to-end las hace Valentino en Slicer con casos reales y pasa la
salida de la consola / del panel "Detalle técnico". Sí se puede verificar que
los archivos compilen (`python -m py_compile`) y pyflakes.
