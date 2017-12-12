/****************************************************************************
** Meta object code from reading C++ file 'decorate_base.h'
**
** Created by: The Qt Meta Object Compiler version 63 (Qt 4.8.5)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../../decorate_base.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'decorate_base.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 63
#error "This file was generated using the moc from 4.8.5. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_ExtraMeshDecoratePlugin[] = {

 // content:
       6,       // revision
       0,       // classname
       0,    0, // classinfo
       2,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       1,       // signalCount

 // signals: signature, parameters, type, tag, flags
      25,   24,   24,   24, 0x05,

 // slots: signature, parameters, type, tag, flags
      57,   48,   24,   24, 0x0a,

       0        // eod
};

static const char qt_meta_stringdata_ExtraMeshDecoratePlugin[] = {
    "ExtraMeshDecoratePlugin\0\0"
    "askViewerShot(QString)\0name,val\0"
    "setValue(QString,vcg::Shotf)\0"
};

void ExtraMeshDecoratePlugin::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        Q_ASSERT(staticMetaObject.cast(_o));
        ExtraMeshDecoratePlugin *_t = static_cast<ExtraMeshDecoratePlugin *>(_o);
        switch (_id) {
        case 0: _t->askViewerShot((*reinterpret_cast< QString(*)>(_a[1]))); break;
        case 1: _t->setValue((*reinterpret_cast< QString(*)>(_a[1])),(*reinterpret_cast< vcg::Shotf(*)>(_a[2]))); break;
        default: ;
        }
    }
}

const QMetaObjectExtraData ExtraMeshDecoratePlugin::staticMetaObjectExtraData = {
    0,  qt_static_metacall 
};

const QMetaObject ExtraMeshDecoratePlugin::staticMetaObject = {
    { &QObject::staticMetaObject, qt_meta_stringdata_ExtraMeshDecoratePlugin,
      qt_meta_data_ExtraMeshDecoratePlugin, &staticMetaObjectExtraData }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &ExtraMeshDecoratePlugin::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *ExtraMeshDecoratePlugin::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *ExtraMeshDecoratePlugin::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_ExtraMeshDecoratePlugin))
        return static_cast<void*>(const_cast< ExtraMeshDecoratePlugin*>(this));
    if (!strcmp(_clname, "MeshDecorateInterface"))
        return static_cast< MeshDecorateInterface*>(const_cast< ExtraMeshDecoratePlugin*>(this));
    if (!strcmp(_clname, MESH_DECORATE_INTERFACE_IID))
        return static_cast< MeshDecorateInterface*>(const_cast< ExtraMeshDecoratePlugin*>(this));
    return QObject::qt_metacast(_clname);
}

int ExtraMeshDecoratePlugin::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QObject::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 2)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 2;
    }
    return _id;
}

// SIGNAL 0
void ExtraMeshDecoratePlugin::askViewerShot(QString _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}
QT_END_MOC_NAMESPACE
