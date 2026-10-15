/* This file is part of MoDELib, distributed under GNU GPL v2. */
#ifndef model_WrappingTabWidget_H_
#define model_WrappingTabWidget_H_
#include <QButtonGroup>
#include <QPushButton>
#include <QResizeEvent>
#include <QStackedWidget>
#include <QScrollArea>
#include <QVBoxLayout>
#include <algorithm>
#include <vector>

namespace model
{
    // Qt's standard tab bar has one row. Keep complete actor names visible
    // in a height-for-width selector above the same persistent control pages.
    class WrappingTabWidget : public QWidget
    {
        class TabBar : public QWidget
        {
            std::vector<QPushButton *> buttons;
            static constexpr int spacing = 4;

            int arrange(int availableWidth, bool place) const
            {
                int x = 0;
                int y = 0;
                int rowHeight = 0;
                for(auto *button : buttons)
                {
                    const QSize size = button->sizeHint();
                    if(x > 0 && x + size.width() > availableWidth)
                    {
                        x = 0;
                        y += rowHeight + spacing;
                        rowHeight = 0;
                    }
                    if(place)
                    {
                        button->setGeometry(x, y, size.width(), size.height());
                    }
                    x += size.width() + spacing;
                    rowHeight = std::max(rowHeight, size.height());
                }
                return y + rowHeight;
            }

          protected:
            void resizeEvent(QResizeEvent *event) override
            {
                QWidget::resizeEvent(event);
                arrange(width(), true);
            }

          public:
            QButtonGroup *const group;

            explicit TabBar(QWidget *parent) : QWidget(parent), group(new QButtonGroup(this))
            {
                QSizePolicy policy(QSizePolicy::Expanding, QSizePolicy::Preferred);
                policy.setHeightForWidth(true);
                setSizePolicy(policy);
                setAccessibleName("Actor tabs");
                group->setExclusive(true);
            }

            void addTab(const QString &name)
            {
                auto *button = new QPushButton(name, this);
                button->setCheckable(true);
                button->setAutoDefault(false);
                group->addButton(button, static_cast<int>(buttons.size()));
                buttons.push_back(button);
                button->setChecked(buttons.size() == 1);
                button->show();
                arrange(width(), true);
                updateGeometry();
            }

            int heightForWidth(int width) const override
            {
                return arrange(width, false);
            }

            QSize minimumSizeHint() const override
            {
                QSize size;
                for(const auto *button : buttons)
                {
                    size = size.expandedTo(button->sizeHint());
                }
                return size;
            }

            QSize sizeHint() const override
            {
                const int preferredWidth = std::max(320, minimumSizeHint().width());
                return QSize(preferredWidth, heightForWidth(preferredWidth));
            }
        };

        TabBar *const bar;
        QStackedWidget *const pages;

      public:
        explicit WrappingTabWidget(QWidget *parent = nullptr)
            : QWidget(parent), bar(new TabBar(this)), pages(new QStackedWidget(this))
        {
            auto *layout = new QVBoxLayout(this);
            layout->setContentsMargins(0, 0, 0, 0);
            layout->setSpacing(4);
            layout->addWidget(bar);
            layout->addWidget(pages, 1);
            connect(bar->group, &QButtonGroup::idClicked, pages, &QStackedWidget::setCurrentIndex);
        }

        int addTab(QWidget *page, const QString &label)
        {
            auto *scroll = new QScrollArea(this);
            scroll->setWidgetResizable(true);
            scroll->setFrameShape(QFrame::NoFrame);
            scroll->setWidget(page);
            const int index = pages->addWidget(scroll);
            bar->addTab(label);
            return index;
        }

        int currentIndex() const
        {
            return pages->currentIndex();
        }
    };
} // namespace model
#endif
