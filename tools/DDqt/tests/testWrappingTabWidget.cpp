#include <WrappingTabWidget.h>
#include <DDqtControls.h>
#include <QLineEdit>
#include <QApplication>
#include <QScrollArea>
#include <QScrollBar>
#include <QSpinBox>
#include <QComboBox>
#include <QSlider>
#include <QWheelEvent>
#include <set>
#include <stdexcept>
#include <iostream>

int main(int argc, char **argv)
{
    QApplication app(argc, argv);
    try
    {
        model::controls::WheelGuard guard;
        app.installEventFilter(&guard);
        QScrollArea scroll;
        auto* body=new QWidget;
        body->setMinimumSize(200,2000);
        auto* spin=new QSpinBox(body);
        spin->setValue(10);
        auto* combo=new QComboBox(body);
        combo->addItems({"one","two","three"});
        combo->setCurrentIndex(1);
        auto* slider=new QSlider(Qt::Horizontal,body);
        slider->setValue(30);
        scroll.setWidget(body);
        scroll.resize(300,300);
        scroll.show();
        app.processEvents();
        for(QWidget* control:std::vector<QWidget*>{spin,spin->findChild<QLineEdit*>(),combo,slider})
        {
            for(bool focused:{false,true})
            {
                scroll.verticalScrollBar()->setValue(200);
                if(focused)
                {
                    control->setFocus();
                }
                else
                {
                    control->clearFocus();
                }
                const auto global=control->mapToGlobal(QPoint(5,5));
                QWheelEvent event(QPointF(5,5),global,QPoint(),QPoint(0,-120),Qt::NoButton,Qt::NoModifier,Qt::NoScrollPhase,false);
                QApplication::sendEvent(control,&event);
                if(spin->value()!=10 || combo->currentIndex()!=1 || slider->value()!=30
                   || scroll.verticalScrollBar()->value()<=200)
                {
                    throw std::runtime_error("Wheel edited a control or failed to scroll the panel");
                }
            }
        }
        spin->stepUp();
        if(spin->value()!=11)
        {
            throw std::runtime_error("Explicit numeric editing was disabled");
        }
        QWidget panel;
        QCheckBox visibility(&panel);
        auto *group = model::controls::box("Test", &panel, "Group help", &visibility);
        auto *edit = new QLineEdit();
        model::controls::field(model::controls::grid(group), 0, "Value", edit, "Detailed field help");
        if(group->isChecked() || edit->isEnabled() || !visibility.isHidden())
        {
            throw std::runtime_error("Initial group state is incorrect");
        }
        group->setChecked(true);
        if(!visibility.isChecked() || !edit->isEnabled())
        {
            throw std::runtime_error("Group checkbox did not enable its controls");
        }
        visibility.setChecked(false);
        if(group->isChecked() || edit->isEnabled())
        {
            throw std::runtime_error("Programmatic visibility did not synchronize the group");
        }
        const auto *label = group->findChild<QLabel *>();
        if(!label || label->buddy() != edit || edit->toolTip().isEmpty() || label->toolTip().isEmpty())
        {
            throw std::runtime_error("Missing field label or tooltip");
        }
        auto *panelLayout = new QVBoxLayout(&panel);
        panelLayout->addWidget(group);
        auto *nested = new model::controls::GroupBox("Nested", group);
        nested->setCheckable(true);
        nested->setChecked(false);
        auto *nestedLayout = new QVBoxLayout(nested);
        auto *nestedEdit = new QLineEdit("kept value");
        nestedLayout->addWidget(nestedEdit);
        model::controls::grid(group)->addWidget(nested, 1, 0, 1, 2);
        auto *intentionallyHidden = new QLabel("Hidden option", group);
        model::controls::grid(group)->addWidget(intentionallyHidden, 2, 0);
        intentionallyHidden->hide();
        panel.show();
        app.processEvents();
        const int collapsedHeight = group->sizeHint().height();
        if(edit->isVisible())
        {
            throw std::runtime_error("Unchecked group contents are visible");
        }
        group->setChecked(true);
        app.processEvents();
        if(!edit->isVisible() || group->sizeHint().height() <= collapsedHeight || nestedEdit->isVisible())
        {
            throw std::runtime_error("Expansion failed or expanded an unchecked nested group");
        }
        nested->setChecked(true);
        group->setChecked(false);
        app.processEvents();
        if(edit->isVisible() || nestedEdit->isVisible() || group->sizeHint().height() != collapsedHeight)
        {
            throw std::runtime_error("Collapse did not reclaim space");
        }
        group->setChecked(true);
        app.processEvents();
        if(!nestedEdit->isVisible() || nestedEdit->text() != "kept value" || !intentionallyHidden->isHidden())
        {
            throw std::runtime_error("Collapse lost nested state, values or intentional visibility");
        }
        model::WrappingTabWidget tabs;
        for(const auto *name : {"Nodes", "Segments", "Loops", "Inclusions", "GlidePlanes", "Quadrature", "Chart",
                                "Fields", "SurfaceMesh", "SlipSystems"})
        {
            tabs.addTab(new QWidget(), name);
        }
        tabs.resize(340, 500);
        tabs.show();
        app.processEvents();
        auto buttons = tabs.findChildren<QPushButton *>();
        const auto checkLayout = [&]()
        {
            std::set<int> rows;
            for(const auto *button : buttons)
            {
                rows.insert(button->y());
                if(button->width() < button->sizeHint().width() ||
                   !button->parentWidget()->rect().contains(button->geometry()))
                {
                    throw std::runtime_error("Actor tab label is clipped");
                }
                for(const auto *other : buttons)
                {
                    if(button != other && button->geometry().intersects(other->geometry()))
                    {
                        throw std::runtime_error("Actor tab buttons overlap");
                    }
                }
            }
            return rows.size();
        };
        const auto narrowRows = checkLayout();
        if(narrowRows < 2)
        {
            throw std::runtime_error("Actor tabs did not wrap");
        }
        buttons[8]->click();
        if(tabs.currentIndex() != 8 || !buttons[8]->isChecked() || buttons[0]->isChecked())
        {
            throw std::runtime_error("Actor tab selection failed");
        }
        tabs.resize(1000, 500);
        app.processEvents();
        if(checkLayout() >= narrowRows || tabs.currentIndex() != 8)
        {
            throw std::runtime_error("Resize lost selection or failed to reflow");
        }
        tabs.resize(240, 500);
        app.processEvents();
        checkLayout();
        std::cout << "PASS: actor tab wrapping, full labels, selection and resize\n";
        return 0;
    }
    catch(const std::exception &error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
