Introducing External Field
==========================

The LFRic infrastructure provides field objects which hold field data in a form
useful to the infrastucture. Other models, libraries and tooling will hold
field data in a different form, useful to their needs. Sometimes it is
necessary to transfer data between these incompatible representations when
linking models together.

The mechanism on offer to handle these occasions is the "External field."

Design
~~~~~~

The infrastructure goes to some lengths to protect the raw field data from
interference. This provides a lot of value in the form of reliability and
maintainability.

Transforming field data for an external user must, by definition, compromise
that protection. The concept of external fields exists to minimise the scope of
that vulnerability.

Each external field is, in fact, a mapping between an LFRic field and a
non-LFRic field. It contains any encapsulation breakage within itself, thus
preventing the internal data of fields from being exposed.

.. figure:: images/external_field_class.svg

    Class diagram showing an external field linking an LFRic field with
    another type of field.

Although the class diagram above shows a mapping between two field classes
there is no reason why the exeternal field need be an actual object. In the
case of coupling through OASIS the external "field" is, in fact, a call to the
transmit or receive functions of OASIS.

Implementing an External Field
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Developing a new external field should not be too hard. The most complicated
part is working out how the external field data is presented.

Start by deriving a new class from the abstract external field class.

Due to the way Fortran handles object construction plain initialiser routines
are used instead but they can be treated like constructors.

The initialiser of the new class should take an LFRic field pointer and some
unique reference to the external target. This may be the target's field
representation (class or primitive array) or some address or ID string.
Whatever is appropriate.

The initialiser should call the abstract initialiser with the LFRic field
pointer. The external reference should be stored for future use.

The two deferred methods (`copy_from_lfric` and `copy_to_lfric`) then need to
be implemented. This is where the work happens.

They both use the field proxy object to gain access to the field's contents and
either populate the external field with that content or update the content from
the external field.
