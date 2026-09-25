#-----------------------------------------------------------------------------
set(MODULE_NAME CranioPlan)

#-----------------------------------------------------------------------------
# Los archivos de CranioPlanLib/ TIENEN que estar listados aca.
#
# Sin esto, al empaquetar CranioPlan como extension de Slicer se copia solo
# CranioPlan.py, la extension se instala sin la libreria, y el modulo revienta
# al arrancar con "ModuleNotFoundError: CranioPlanLib". Andando por
# "Additional module paths" no se nota, porque ahi se usa la carpeta tal cual
# esta en el disco; el problema aparece recien al distribuirlo.
set(MODULE_PYTHON_SCRIPTS
  ${MODULE_NAME}.py
  CranioPlanLib/__init__.py
  CranioPlanLib/Comun.py
  CranioPlanLib/BloqueC.py
  CranioPlanLib/BloqueA.py
  CranioPlanLib/BloqueF.py
  CranioPlanLib/BloqueG.py
  CranioPlanLib/PuenteFG.py
  )

set(MODULE_PYTHON_RESOURCES
  Resources/Icons/${MODULE_NAME}.png
  Resources/UI/${MODULE_NAME}.ui
  )

#-----------------------------------------------------------------------------
slicerMacroBuildScriptedModule(
  NAME ${MODULE_NAME}
  SCRIPTS ${MODULE_PYTHON_SCRIPTS}
  RESOURCES ${MODULE_PYTHON_RESOURCES}
  WITH_GENERIC_TESTS
  )

#-----------------------------------------------------------------------------
if(BUILD_TESTING)

  # Register the unittest subclass in the main script as a ctest.
  # Note that the test will also be available at runtime.
  slicer_add_python_unittest(SCRIPT ${MODULE_NAME}.py)

  # Additional build-time testing
  add_subdirectory(Testing)
endif()
