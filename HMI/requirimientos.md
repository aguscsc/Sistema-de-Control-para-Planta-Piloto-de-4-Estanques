# Requerimientos de HMI

## Pantallas:
    
    - *Pantalla de Inicio*: Debe existir una pantalla de inico con usuario y contraseña (gestión de permisos)
    - *P&ID*: debe tener el diagrama P&ID de la planta, los gráficos de tendencia de los procesos de interés (Nivel de estanques, flujo)
    - *Pantalla de sintonización*: Debe existir una pantalla que permita sintonizar los lazos existentes

## Elementos Fijos
    
    - *Parada de Emergencia*: En colores saturados
    - *Alarmas*: Alarmas (LL, L, H, HH)
    - *Menú*: cotiene pantallas disponibles
    - *Barra de inofrmación*: Contiene alarmas, hora, fecha y usuario

## Registros
        
    - *Todas las alarmas deben ser registradas*

## Permisos
    
    - *Administrador*: Permiso a todas las operaciones.
    - *Operador*: Puesta marcha, Parada de Emergencia,cambiar límites de variables, cambio automático/manual (NO puede cambiar parámetros del controlador).
    _ *Visitante*: Sólo puede observar los procesos

## Aspectos visuales
    
    - *Color de fondo*: Fondo gris claro
    - *Colores saturados*: Reservados sólo para alarmas y para de Emergencia
    - *Animaciones*: Evitar animaciones o movimientos en la pantalla que pudieran distraer al operador
    - *POP-UPS*: No pueden ser escalados
    - *Evitar gradiente*
    - *Diferenciación de estado de botones*: Gris oscuro para selección y no seleccionado gris claro
    
### Textos y Presentación de datos
    
    - *Fuentes y formato*: Las fuentes para los textos deben ser Arial, con tamaño de fuentes 10 y 12 dependiendo de la densidad de la pantalla. El atributo ennegrecido se utiliza dependiendo de la densidad y cuanto se requiera destacar un mensaje.
    - 
    - *Barra analógica*: Mostrar barra de comparación entre Sp, Pv y Cv

